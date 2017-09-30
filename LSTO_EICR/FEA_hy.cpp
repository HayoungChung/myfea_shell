#include "FEA_hy.h"
#include <fstream>
#include <cmath>

FEAMesh::FEAMesh(const double* _Lxy, const int* _exy, bool isOverlaid){
    Lxy = _Lxy; exy = _exy;
    this->isOverlaid = isOverlaid;
    get_Mesh(3);
    // initialize properties
    areafraction.reserve(nELEM);
    centroid = MatrixXd::Zero(nELEM, 3);
    ElemArea.reserve(nELEM);
    set_ElemArea();
    ComputeCentroids();
}

void FEAMesh::set_ElemArea(){
    double multi;
    VectorXi elem_id;
    multi = (isOverlaid == true) ? 1.0 : 0.5;
    for (unsigned int ee = 0; ee < nELEM; ++ee){
        // TODO: assume structured mesh
        elem_id = ELEM.row(ee);
        VectorXd l1 = NODE.row(elem_id(1)) - NODE.row(elem_id(0));
        VectorXd l2 = NODE.row(elem_id(2)) - NODE.row(elem_id(0));
        ElemArea[ee] = std::abs(l1(0)*l2(1) - l2(0)*l1(1))*multi; 
        // std::cout << ee << " : " << ElemArea[ee] << std::endl;
    }
}

void FEAMesh::get_Mesh(int nVertices){
    // static allocation is prohibited unless the size is known. no variables accepted
    double lx = Lxy[0], ly = Lxy[1];
    int ex = exy[0], ey = exy[1];
    int nxQ = ex+1, nyQ = ey + 1; // Q4 mesh

    unsigned int nNODEQ = nxQ*nyQ; 
    unsigned int nELEMQ = ex*ey;
    bool & isOverlaid = this->isOverlaid;
    // Meshgrid (following M2DO convention)
    MatrixXd ndx = VectorXd::LinSpaced(nxQ,0,lx).replicate(1,nyQ);
    MatrixXd ndy = RowVectorXd::LinSpaced(nyQ,0,ly).replicate(nxQ,1);

    // Q4 meshes (can be used when overrided mesh)
    MatrixXd NODEQ; NODEQ.setZero(nNODEQ, 3);
    NODEQ.col(0) = Map<MatrixXd>(ndx.data(),nNODEQ,1);
    NODEQ.col(1) = Map<MatrixXd>(ndy.data(),nNODEQ,1);

    // Matrix<double, nELEMQ, npe+1> ELEMQ; ELEMQ.fill(0);
    MatrixXi ELEMQ; ELEMQ.setZero(nELEMQ, npe+1);
    for (int ii = 0; ii < ex; ++ii){
        for (int jj = 0; jj < ey; ++jj){
            int eid = jj*ex + ii;
            int nid = jj*nxQ + ii;
            ELEMQ(eid,0) = nid;
            ELEMQ(eid,1) = nid + 1;
            ELEMQ(eid,2) = nid + 1 + nxQ;
            ELEMQ(eid,3) = nid + nxQ;            
        }
    }
    if (nVertices == 4){
        this->ELEM = ELEMQ;
        this->NODE = NODEQ;
        dpn = 2, npe = 4, dpe = 8; // Quadmesh
        this->nELEM = nELEMQ;
        this->nNODE = nNODEQ;
        return;
    }
    if (isOverlaid == true){
        // overlaid
        this->ELEM = ELEMQ;
        this->NODE = NODEQ;
        dpn = 6; npe = 3; dpe = 18; 
        this->nELEM = nELEMQ;
        this->nNODE = nNODEQ;
        return;
    }

    // Triangular mesh
    nELEM = ex*ey*4; // Lattice-shaped mesh (no override)
    // Matrix<double, nELEM, npe> ELEM; ELEM.fill(0);
    MatrixXi ELEM; ELEM.setZero(nELEM,npe);

    nNODE = nNODEQ + nELEMQ;
    // Matrix<double, nNODE, 3> NODE; 
    MatrixXd NODE; NODE.setZero(nNODE,3);

    NODE.topRows(nNODEQ) = NODEQ;

    // change the order to get triangle
    Matrix<int, 4, 3> Tri_order; 
    Tri_order << 4, 1, 5, 1, 2, 5, 2, 3, 5, 3, 4, 5;
    Tri_order.array() -= 1;

    for (unsigned int ee = 0; ee < nELEMQ; ++ ee){
        int nid = nNODEQ + ee;
        VectorXi elem_id = ELEMQ.row(ee);
        // compute centeroid
        double x_sum =0, y_sum = 0, z_sum = 0;
        
        for (unsigned int cc = 0; cc < 4; ++cc){
            x_sum += NODEQ(elem_id(cc),0);
            y_sum += NODEQ(elem_id(cc),1);
            z_sum += NODEQ(elem_id(cc),2);
        }
        NODE(nNODEQ+ee,0) = x_sum/4;
        NODE(nNODEQ+ee,1) = y_sum/4;
        NODE(nNODEQ+ee,2) = z_sum/4;

        // change the order to get triangle
        for (int gg = 0; gg < 4; ++gg){
            int eid_T = 4*ee+gg;
            ELEM(eid_T,0) = elem_id(Tri_order(gg,0));

            ELEM(eid_T,1) = elem_id(Tri_order(gg,1));
            // ELEM(eid_T,2) = ELEMQ(ee,Tri_order(gg,2))
            ELEM(eid_T,2) = nNODEQ+ee;
        }        
    }
    this->ELEM = ELEM;
    this->NODE = NODE;
    dpn = 6, npe = 3, dpe = 18; // tri shell

}

std::vector<int> FEAMesh::get_nodeID(double posX, double posY, double Xtol, double Ytol){
    Array<bool, Dynamic, 1> Xid, Yid;
    Xid = (abs(NODE.col(0).array() - posX) < Xtol);
    Yid = (abs(NODE.col(1).array() - posY) < Ytol);
    std::vector<int> NodeId;
    for (unsigned int nn = 0; nn < NODE.rows(); ++nn){
        if (Xid(nn) == false || Yid(nn) == false) continue;
        NodeId.push_back(nn);
    }
    return NodeId;
}

std::vector<int> FEAMesh::get_dof(int direction, std::vector<int> & ID_nodes){
    std::vector<int> Dofs;
    // Dofs.reserve(ID_nodes.size());
    for (unsigned int ii = 0; ii < ID_nodes.size(); ++ii){
        Dofs.push_back( ID_nodes[ii]*dpn + direction ); 
    }     
    return Dofs;
}

void FEAMesh::get_BCid(std::vector<int> Dofs){
    std::vector<int> & bc_checked = Dofs;
    
    std::sort (bc_checked.begin(), bc_checked.end());
    std::vector<int>::iterator last = std::unique (bc_checked.begin(), bc_checked.end());
    bc_checked.erase(last, bc_checked.end());

    this->BCid.setZero(bc_checked.size());
    for (unsigned int nn = 0; nn < bc_checked.size(); ++nn){
        this->BCid(nn) = bc_checked[nn];
    }
}
 
void FEAMesh::set_Force(int direction, std::vector<int>& ID_nodes, double val, MatrixXd & force_fix){
    // if (ID_nodes.size() == 0) return;
    for (unsigned int ii = 0; ii < ID_nodes.size(); ++ii){
        force_fix(ID_nodes[ii], direction) += val;
    }
}

void FEAMesh::set_Force(std::vector<Material_ABD> & material, MatrixXd eps0, MatrixXd kappa0, MatrixXd & force_NM){
    VectorXd NM = VectorXd::Zero(dpn); 
    for (unsigned int ee = 0; ee < ELEM.rows(); ++ee ){
        Vector3d eps_ = eps0.row(ee), kappa_ = kappa0.row(ee);
        NM << material[ee].Amat*eps_, material[ee].Dmat*kappa_;
        force_NM.row(ee) = NM.transpose();
    }
}

void FEAMesh::to_vtk(MatrixXd & u, const char* str){
    int cell_type;
    switch (npe){
        case 3: //elemType = 'triangle';
                cell_type = 5;
        case 4: // elemType = 'quadrilateral';
                cell_type = 9;
    }
    int npe_draw = npe;
    if (isOverlaid == true) {
        cell_type = 9;
        npe_draw = 4;
    }
    int nNODE = NODE.rows();
    int nELEM = ELEM.rows();

    std::ofstream myfile;
    myfile.open(str);
    myfile << "# vtk DataFile Version 3.0\n" ;
	myfile << "vtk output\n" ;
	myfile << "ASCII\n\n" ;
	myfile << "DATASET UNSTRUCTURED_GRID" ;
	myfile << "\nPOINTS " << nNODE << " double" ;

    // NODE
    myfile << NODE;

    // ELEMENT
    int nCellPoints = nELEM * (npe_draw + 1);
	myfile << "\n\nCELLS " << nELEM << " " << nCellPoints;
	
	for (int i = 0 ; i < nELEM ; ++i) {
        myfile << "\n" << npe_draw << " " << ELEM.row(i) << " ";
    }

	myfile << "\n\nCELL_TYPES " << nELEM ;
	// int cell_type = (spacedim == 2) ? 9 : 12 ; (12: 3d quad)
	
	for (int i = 0 ; i < nELEM ; ++i) {
		myfile << "\n" << cell_type ;
	}

    if (u.rows() == nNODE){
        myfile << "\n\nPOINT_DATA " << nNODE ;
        if (u.cols() == 1){ //scalar
            myfile << "\n\nSCALARS u double \n" ;
            for (int i = 0 ; i < nNODE ; ++i) {
                myfile << u(i) << " \n";
            }
        }
        else {
            myfile << "\n\nVECTORS u double \n" ;
            for (int i = 0 ; i < nNODE ; ++i) {
                myfile << u(i,0) << " " << u(i,1) << " " << u(i,2) << " \n";
            }
        }
    }
    else if (u.rows() == nELEM){
        myfile << "\n\nCELL_DATA " << nELEM ;
        if (u.cols() == 1){ //scalar
            myfile << "\n\nSCALARS u double \n" ;
            for (int i = 0 ; i < nELEM ; ++i) {
                myfile << u(i,0) << " \n";
            }
        }
        else {
            myfile << "\n\nVECTORS u double \n" ;
            for (int i = 0 ; i < nELEM ; ++i) {
                myfile << u(i,0) << " " << u(i,1) << " " << u(i,2) << " \n";
            }
        }
    }

	myfile.close();
}

void FEAMesh::to_vtk(MatrixXd & u){
    this->to_vtk(u, "mesh_without_field.vtk");
}

void FEAMesh::to_vtk(){
    int cell_type;
    switch (npe){
        case 3: //elemType = 'triangle';
                cell_type = 5;
        case 4: // elemType = 'quadrilateral';
                cell_type = 9;
    }
    int npe_draw = npe;
    if (isOverlaid == true) {
        cell_type = 9;
        npe_draw = 4;
    }
    int nNODE = NODE.rows();
    int nELEM = ELEM.rows();

    std::ofstream myfile;
    myfile.open("mesh_without_fields.vtk");
    myfile << "# vtk DataFile Version 3.0\n" ;
	myfile << "vtk output\n" ;
	myfile << "ASCII\n\n" ;
	myfile << "DATASET UNSTRUCTURED_GRID" ;
	myfile << "\nPOINTS " << nNODE << " double" ;

    // NODE
    myfile << NODE;

    // ELEMENT
    int nCellPoints = nELEM * (npe_draw + 1);
	myfile << "\n\nCELLS " << nELEM << " " << nCellPoints;
	
	for (int i = 0 ; i < nELEM ; ++i) {
        myfile << "\n" << npe_draw << " " << ELEM.row(i) << " ";
    }

	myfile.close();
}

void FEAMesh::ComputeCentroids(){
    // Overlaid case 
    MatrixXi elem_id;
    // double multiple = (isOverlaid == true) ? 4.0 : 3.0;
    int nE = ELEM.cols();
    for (unsigned int ee = 0; ee < nELEM; ++ee){
           elem_id = ELEM.row(ee);
           
           for (unsigned int mmm = 0; mmm < nE; ++mmm){
               centroid(ee,0) += NODE(elem_id(mmm),0)/double(nE);
               centroid(ee,1) += NODE(elem_id(mmm),1)/double(nE);
               centroid(ee,2) += NODE(elem_id(mmm),2)/double(nE);
           }
    }    
}