

/****** Functions related to state equation of unknown "u" ******/

double func_a(double x, double y);
double func_c(double x, double y);
double func_d(double x, double y);
double func_g(double x, double y);
double uD(double x, double y);






/******* Functions related to adjoint equation of unknown "xi" *********/

double func_a_adjoint(double x, double y);
double func_c_adjoint(double x, double y);
double func_d_adjoint(double x, double y);
double xiD(double x, double y);



/******* Function to compute L2 projection of normal derivative on interface for a given solution ************/
int computeInterfaceNormalDerivative(	const my_mesh::MeshData &mesh, 
					const arma::vec &u,
					arma::vec &du_dnu_1, arma::vec &du_dnu_2);
