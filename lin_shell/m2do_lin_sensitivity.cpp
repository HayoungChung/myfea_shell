#include "./m2do_lin_sensitivity.h"

// Constructor
SensitivityAnalysis::SensitivityAnalysis(FEAMesh& feaMesh_, std::vector<GptsCompl> &gptsCompls_) : gptsCompls(gptsCompls_), feaMesh(feaMesh_)
{
    int sens_order;
    if (feaMesh.isOverlaid == true)
    {
        sens_order = 2;
    }
    else
    {
        sens_order = 1;
    }
    // Declaring variables
    int nel = feaMesh.ELEM.rows();   // Total number of elements
    int ngauss = pow(sens_order, 2); // Total number of gauss points

    // Resizing sensitivities in function of nel and ngauss
    sensitivities.resize(nel);
    for (int i = 0; i < nel; i++)
    {
        sensitivities[i].dx.resize(ngauss);
        sensitivities[i].dx1.resize(ngauss);
        sensitivities[i].dx2.resize(ngauss);
        sensitivities[i].dxcoord.resize(ngauss);
        for (int j = 0; j < ngauss; j++)
        {
            sensitivities[i].dxcoord[j].resize(spacedim);
        }

        // Initialising list of excluded elements.
        sensitivities[i].isExcluded = false;
    }
    GptsCompl2sensitivities();
}

void SensitivityAnalysis::GptsCompl2sensitivities()
{
    int sens_order;
    if (feaMesh.isOverlaid == true)
    {
        sens_order = 2;
    }
    else
    {
        sens_order = 1;
    }
    // Declaring variables
    int nel = feaMesh.ELEM.rows();   // Total number of elements
    int ngauss = pow(sens_order, 2); // Total number of gauss points

    for (int ee = 0; ee < nel; ee++)
    {
        for (int gg = 0; gg < ngauss; gg++)
        {
            sensitivities[ee].dx[gg] = gptsCompls[ee * ngauss + gg].sens;
            sensitivities[ee].dxcoord[gg][0] = gptsCompls[ee * ngauss + gg].x;
            sensitivities[ee].dxcoord[gg][1] = gptsCompls[ee * ngauss + gg].y;
        }
    }
}

// LEAST SQUARES FUNCTIONALITIES

// Functionality to compute sensitivities at the level-set boundaries.
// template <int spacedim, int dim, int order, class Mesh, class Physics, class Study, class LSM, int sens_order>
void SensitivityAnalysis::ComputeBoundarySensitivities(double radius, std::vector<double> bPoint, double areamin)
{
    // DECLARING VARIABLES
    int sens_order;
    if (feaMesh.isOverlaid == true)
    {
        sens_order = 2;
    }
    else
    {
        sens_order = 1;
    }
    // Scalars and vectors
    int nel = feaMesh.ELEM.rows();          // Total number of elements
    int ngauss = pow(sens_order, spacedim); // Total number of gauss points
    double rSqd = radius * radius;          // Squared radius.
    double dist;                            // Squared distance from the Gauss point to the element i.
    // std::vector<double> signedDistGP(ngauss,0); // Signed distance function values at the Gauss points.

    // Least square information object.
    LeastSquares leastsq;

    // FUNCTION BODY

    // Finding Gauss points inside the boundary and the subdomain defined by radius
    for (int i = 0; i < nel; i++)
    {
        // Looking at nearby elements
        int aux_out = 0;
        for (int j = 0; j < spacedim; j++)
        {
            // If the element is too far from boundary point -> aux_out = 1.
            if ((feaMesh.centroid(i, j) > (bPoint[j] + 1.5 * radius)) || (feaMesh.centroid(i, j) < (bPoint[j] - 1.5 * radius)))
            {
                aux_out = 1;
            }
        }
        // If element is outside the subdomain defined by 1.5radius
        if (aux_out == 1)
        {
            // Element is too far from boundary point
        }
        else // Else element is close enough to be considered
        {
            // Check if element is inside the boundary
            if (feaMesh.areafraction[i] > areamin)
            {
                // // Checking if element is a cut element and then check signed distance function at the Gauss points.
                // if (feaMesh.areafraction[i] < 1)
                // {
                //     signedDistGP = ComputeSignedDistanceAtGaussPoints(mesh, i);
                // }

                // For each Gauss point of the ith element
                for (int j = 0; j < ngauss; j++)
                {
                    // // Take into account only Gauss points with positivite signed distance.
                    // if (signedDistGP[j] >= 0)
                    // {
                    // Compute squared distance between current Gauss point and boundary point
                    dist = 0.0;
                    for (int k = 0; k < spacedim; k++)
                    {
                        // Squared distance in each space dimension
                        dist += pow(bPoint[k] - sensitivities[i].dxcoord[j][k], 2);
                    }

                    // If squared distance is less than squared radius, save information
                    if (dist < rSqd)
                    {
                        // Grouping information in object
                        leastsq.distGauss = sqrt(dist);                      // Save distance
                        leastsq.areaFractionGauss = feaMesh.areafraction[i]; // Save area fraction
                        leastsq.elementNumber = i;                           // Save element number
                        leastsq.gaussPointNumber = j;                        // Save Gauss point number
                        leastsq.coord = sensitivities[i].dxcoord[j];         // Coordinates
                        // Storing information in class
                        leastsquares.push_back(leastsq);
                    }
                    // }
                }
            }
        }
    }

    // If not enough Gauss points inside subdomain, we are on an island.
    if (leastsquares.size() < 10)
    {
        // Sensitivities at the boundaries for islands are zero.
        boundarySens.push_back(0.0);
        return;
    }

    // SOLVING LEAST SQUARES

    // Solve least squares and obtain sensitivity at the boundary
    double B = SolveLeastSquares(leastsquares, bPoint);

    // Store sensitivity
    boundarySens.push_back(B);

    // Clear leastsquares vector.
    leastsquares.clear();
}

// Functionality to solve least squares problem for 2D and 3D cases.
double SensitivityAnalysis::SolveLeastSquares(std::vector<LeastSquares>& leastsquares, std::vector<double>& bPoint)
{

    // Number of Gauss points inside subdomain defined by radius.
    int npoints = leastsquares.size();

    // Size of elements in least squares' basis function.
    int basis[2] = {6, 10};
    // Selecting size according to dimensionality of the problem.
    int n = basis[spacedim - 2];

    // Declaring variables.
    double A[n * npoints];
    double B[npoints];
    double Bmin = +1e6;
    double Bmax = -1e6;
    double sens;

    // Building least squares problem.
    if (spacedim == 2) // For 2D interpolation.
    {
        for (int i = 0; i < npoints; i++)
        {
            // Weight function by inverse distance.
            double lsweight = leastsquares[i].areaFractionGauss / leastsquares[i].distGauss;

            // Relative x and y coordinates.
            double xb = leastsquares[i].coord[0] - bPoint[0];
            double yb = leastsquares[i].coord[1] - bPoint[1];

            // Storing weighted distances information.
            A[i] = lsweight;
            A[npoints + i] = xb * lsweight;
            A[(2 * npoints) + i] = yb * lsweight;
            A[(3 * npoints) + i] = xb * yb * lsweight;
            A[(4 * npoints) + i] = xb * xb * lsweight;
            A[(5 * npoints) + i] = yb * yb * lsweight;

            // Sensitivity at the current point.
            sens = sensitivities[leastsquares[i].elementNumber].dx[leastsquares[i].gaussPointNumber];

            // Storing weighted sensitivity.
            B[i] = sens * lsweight;

            if (sens > Bmax)
                Bmax = sens;
            else if (sens < Bmin)
                Bmin = sens;
        }
    }
    else if (spacedim == 3) // For 3D interpolation.
    {
        for (int i = 0; i < npoints; i++)
        {
            // Weight function by inverse distance.
            double lsweight = leastsquares[i].areaFractionGauss / leastsquares[i].distGauss;

            // Relative x, y and z coordinates.
            double xb = leastsquares[i].coord[0] - bPoint[0];
            double yb = leastsquares[i].coord[1] - bPoint[1];
            double zb = leastsquares[i].coord[2] - bPoint[2];

            // Storing weighted distances information.
            A[i] = lsweight;
            A[npoints + i] = xb * lsweight;
            A[(2 * npoints) + i] = yb * lsweight;
            A[(3 * npoints) + i] = zb * lsweight;
            A[(4 * npoints) + i] = xb * yb * lsweight;
            A[(5 * npoints) + i] = xb * zb * lsweight;
            A[(6 * npoints) + i] = yb * zb * lsweight;
            A[(7 * npoints) + i] = xb * xb * lsweight;
            A[(8 * npoints) + i] = yb * yb * lsweight;
            A[(9 * npoints) + i] = zb * zb * lsweight;

            // Sensitivity at the current point.
            sens = sensitivities[leastsquares[i].elementNumber].dx[leastsquares[i].gaussPointNumber];

            // Storing weighted sensitivity.
            B[i] = sens * lsweight;

            if (sens > Bmax)
                Bmax = sens;
            else if (sens < Bmin)
                Bmin = sens;
        }
    }

    // SOLVE LEAST SQUARES WITH LAPACK

    // Auxiliary variables.
    int m = npoints, nrhs = 1, lda = npoints, ldb = npoints, info, lwork;
    double wkopt;
    double *work;

    // Executable statements.
    // Query and allocate the optimal workspace.
    lwork = -1;
    dgels_("No transpose", &m, &n, &nrhs, A, &lda, B, &ldb, &wkopt, &lwork, &info);

    // Allocating memory.
    lwork = (int)wkopt;
    work = (double *)malloc(lwork * sizeof(double));

    // Solve equations A*X = B (B is overwritten with the solution X).
    dgels_("No transpose", &m, &n, &nrhs, A, &lda, B, &ldb, work, &lwork, &info);
    free(work);

    double abmax = (Bmax > 0.0) ? 10.0 * Bmax : 0.1 * Bmax;
    double abmin = (Bmin < 0.0) ? 10.0 * Bmin : 0.1 * Bmin;

    // Use weighted average.
    if (B[0] > abmax || B[0] < abmin)
    {
        B[0] = 0.0;
        double wgt = 0.0;

        for (int i = 0; i < npoints; i++)
        {
            // Weight function by inverse distance.
            double lsweight = leastsquares[i].areaFractionGauss / leastsquares[i].distGauss;

            B[0] += sensitivities[leastsquares[i].elementNumber].dx[leastsquares[i].gaussPointNumber] * lsweight;
            wgt += lsweight;
        }
        B[0] /= wgt;
    }
    else if (B[0] > Bmax)
    {
        B[0] = Bmax;
    }
    else if (B[0] < Bmin)
    {
        B[0] = Bmin;
    }

    // Sensitivity at boundary point.
    sens = B[0];

    return sens;
}
