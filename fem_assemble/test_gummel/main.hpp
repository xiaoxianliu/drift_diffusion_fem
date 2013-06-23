

// Functions associated with equation at hand
double psi_func(double x, double y);
double f_func(double x, double y);


// Apply Dirichlet boundary conditions to linear system (coefficient matrix and right-hand side vector)
int applyDirichletBC(	const my_mesh::MeshData &mesh,
			arma::mat &M,
			arma::vec &rhs);


// Plots
int plotSolution(const my_mesh::MeshData& mesh, const arma::vec u, const std::string& filename);
