// this computes a linear sensitivity 


Sensitivity::Sensitivity(FEAMesh& feaMesh_, std::vector<Material_ABD>& material_, MatrixXd& GU_u_): feaMesh(feaMesh_), material(material_), GU_u(GU_u_){
    
}