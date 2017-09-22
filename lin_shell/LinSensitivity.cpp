#include "./LinSensitivity.h"

LinSensitivity::LinSensitivity(FEAMesh &feaMesh_, std::vector<GptsCompl> &gaussComplSens_) : gaussComplSens(gaussComplSens_), feaMesh(feaMesh_)
{
    this->areafraction = feaMesh.areafraction;
    nGpts = feaMesh.isOverlaid == true ? feaMesh.ELEM.rows() * 4 : feaMesh.ELEM.rows();

    Gpts.resize(nGpts, 2);
    GptsSensitivities.resize(nGpts, 1);

    for (unsigned int nn = 0; nn < nGpts; nn++)
    {
        Gpts(nn, 0) = gaussComplSens[nn].x;
        Gpts(nn, 1) = gaussComplSens[nn].y;
        // Gpts(nn, 2) = gaussComplSens[nn].z;

        GptsSensitivities(nn) = gaussComplSens[nn].sens;
    }
}

std::vector<double> LinSensitivity::get_bptsCompl(MatrixXd &BptsCoords, double radius)
{
    unsigned int WeightFlag = 5;

    nBpts = BptsCoords.size();
    BptsSensitivities.resize(nBpts);
    for (unsigned int nn = 0; nn < nBpts; nn++)
    {
        std::vector<double> PointXYZ(3,0.0);
        PointXYZ[0] = BptsCoords(nn, 0);
        PointXYZ[1] = BptsCoords(nn, 1);
        PointXYZ[2] = BptsCoords(nn, 2);
// 
        // BptsSensitivities[nn] = ComputeBoundaryPointSensitivity2D(PointXYZ, radius, WeightFlag, 0.01); // TOFIX: temporarly, only 2D (i.e., flat) config is considered
    }
    return BptsSensitivities;
}

double LinSensitivity::ComputeBoundaryPointSensitivity2D(std::vector<double> &Pointxy, double Radius, unsigned int WeightFlag, double Tolerance)
{
    unsigned int Counter = 0;
    unsigned int p;
    double el_dist;
    int nGptsElem;
    double temp;
    if (feaMesh.isOverlaid == true)
    {
        p = ((PI * Radius * Radius) / feaMesh.ElemArea[0]) * 16;
        nGptsElem = 4;
    }
    else
    {
        p = ((PI * Radius * Radius) / feaMesh.ElemArea[0]);
        nGptsElem = 1;
    }
    p *= 1.25; // for a conservative estimate
    std::vector<double> Distances(p);
    std::vector<unsigned int> ElementIndices(p);
    std::vector<unsigned int> Indices(p);

    double PointSensitivity = 0.0;

    int CntPoints = 0;

    Vector2d el_cood, gg_cood;
    double gg_dist;
    for (unsigned int ee = 0; ee < feaMesh.ELEM.rows(); ++ee)
    {
        // initial domaincounter
        el_cood << feaMesh.Centeroids(ee, 0), feaMesh.Centeroids(ee, 1);

        if (feaMesh.areafraction[ee] > Tolerance)
        {
            el_dist = std::sqrt(std::pow(Pointxy[0] - el_cood[0], 2) + std::pow(Pointxy[1] - el_cood[1], 2));
            if (el_dist < 1.5 * Radius)
            {
                for (int gg = 0; gg < nGptsElem; ++gg)
                {
                    unsigned int gg_indx = nGptsElem * ee + gg; // ee *= nGptsElem 로 인식하는 것 같다.. eigen lib. 문제 Jun20
                    gg_cood = Gpts.row(gg_indx);
                    gg_dist = std::sqrt(std::pow(Pointxy[0] - gg_cood[0], 2) + std::pow(Pointxy[1] - gg_cood[1], 2));
                    if (gg_dist < Radius)
                    {
                        Distances[CntPoints] = gg_dist;
                        ElementIndices[CntPoints] = ee;
                        Indices[CntPoints] = gg;
                        CntPoints += 1;
                    }
                }
            }
        }
    }

    if (CntPoints < 10)
    {
        std::cout << "a very small island is found at (" << Pointxy[0] << ", " << Pointxy[1] << "): npts = " << CntPoints << std::endl;
        PointSensitivity = 0.0;
        return PointSensitivity;
    }

    MatrixXd A(CntPoints, 6);
    VectorXd b_sens(CntPoints);
    double BMax = 1e20, BMin = -1e20;
    std::vector<double> RelativeCoordinate(2);

    // TOFIX: spacedim = 2, nDual = 1
    for (int ii = 0; ii < CntPoints; ++ii)
    {
        p = ElementIndices[ii] * nGptsElem + Indices[ii]; //Indices[ii];
        switch (WeightFlag)
        {
        case 1:
            temp = 1.0; // least squares
            break;
        case 2:
            temp = 1.0 / Distances[ii];
            break;
        case 3:
            temp = feaMesh.areafraction[ii];
            break;
        case 4:
            temp = feaMesh.areafraction[ii] / Distances[ii];
            break;
        case 5:
            temp = std::sqrt(feaMesh.areafraction[ii] / Distances[ii]);
            break;
        default:
            temp = 1.0;
            std::cout << "Weight Flag should lie in [1, 5]. Using Least Squares.\n";
        }
        for (unsigned int jj = 0; jj < 2; ++jj)
        {
            RelativeCoordinate[jj] = Gpts(p, jj) - Pointxy[jj];
        }
        A(ii, 0) = temp;
        A(ii, 1) = RelativeCoordinate[0] * temp;
        A(ii, 2) = RelativeCoordinate[1] * temp;
        A(ii, 3) = RelativeCoordinate[0] * RelativeCoordinate[1] * temp;
        A(ii, 4) = RelativeCoordinate[0] * RelativeCoordinate[0] * temp;
        A(ii, 5) = RelativeCoordinate[1] * RelativeCoordinate[1] * temp;

        b_sens(ii) = -GptsSensitivities(p, 0) * temp; // linear case test
        // if (isLinear == true) b_sens(ii) = -GptsSensitivities (p,0)*temp; // linear case test
        // else b_sens(ii) =  GptsSensitivities (p,3)*temp;
    }
    MatrixXd B_tmp = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b_sens);
    double B = B_tmp(0);
    if ((B > BMax * 10) || (B < BMin * 10))
    {
        B = 0.0;
        temp = 0.0;

        for (int nn = 0; nn < CntPoints; ++nn)
        {
            B += GptsSensitivities(ElementIndices[nn] * nGptsElem + Indices[nn], 3) * feaMesh.areafraction[ElementIndices[nn]];
            temp += feaMesh.areafraction[ElementIndices[nn]];
        }
        PointSensitivity = B / temp;
    }
    else if (B > BMax)
        PointSensitivity = BMax;
    else if (B < BMin)
        PointSensitivity = BMin;
    else
        PointSensitivity = B;

    return PointSensitivity;
}
