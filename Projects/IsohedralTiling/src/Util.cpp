#include <unordered_set>
#include "../include/Util.h"


void loadQuadraticTriangleMeshFromVTKFile(const std::string& filename, Eigen::MatrixXd& V, 
    Eigen::MatrixXi& F, Eigen::MatrixXi& V_quad)
{
    using TV3 = Vector<T, 3>;
    using IV3 = Vector<int, 3>;
    using IV6 = Vector<int, 6>;

    std::ifstream in(filename);
    std::string token;

	while(token != "POINTS")
    {
		in >> token;
        if (in.eof())
            break;
    }
    
    int n_points;
	in >> n_points;

	in >> token; //double

	V.resize(n_points, 3);

	for(int i=0; i<n_points; i++)
        for (int  j = 0; j < 3; j++)
            in >> V(i, j);
        
    while(token != "CELLS")
    {
		in >> token;
        if (in.eof())
            break;
    }
    int n_cells, n_entries;
	in >> n_cells;
	in >> n_entries;
    
    int cell_type;
	std::vector<IV3> faces;
    std::vector<IV6> quad_node_indices;
    // std::cout << n_cells << std::endl;
	for(int i=0; i<n_cells; ++i)
	{
		in >> cell_type;

		if(cell_type == 3)
		{
            IV3 face;
			for(int j = 0; j < 3; j++)
				in >> face[j];
            faces.push_back(face);
		}
        else if (cell_type == 6)
        {
            IV6 quad_nodes;
            for (int j = 0; j < 6; j++)
            {
                in >> quad_nodes[j];
            }
            quad_node_indices.push_back(quad_nodes);
        }
		else
		{
			// type 1 2
			for(int j = 0; j < cell_type; j++)
				in >> token;
		}
	}
    int n_faces = faces.size();
    F.resize(n_faces, 3);
    for (int i = 0; i < n_faces; i++)
    {
        TV3 ei = V.row(faces[i][0]) - V.row(faces[i][1]);
        TV3 ej = V.row(faces[i][2]) - V.row(faces[i][1]);
        if (ej.cross(ei).dot(TV3(0, 0, 1)) < 0)
            F.row(i) = IV3(faces[i][1], faces[i][0], faces[i][2]);
        else
            F.row(i) = faces[i];
    }
    int n_quad_nodes = quad_node_indices.size();
    V_quad.resize(n_quad_nodes, 6);
    for (int i = 0; i < n_quad_nodes; i++)
    {
        V_quad.row(i) = quad_node_indices[i];
    }
    in.close();
}


void loadMeshFromVTKFile3D(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXi& T)
{
    using TV3 = Vector<double, 3>;
    using IV3 = Vector<int, 3>;
    using IV4 = Vector<int, 4>;

    struct IV3Hash
    {
        size_t operator()(const IV3& a) const{
            std::size_t h = 0;
            for (int d = 0; d < 3; ++d) {
                h ^= std::hash<int>{}(a(d)) + 0x9e3779b9 + (h << 6) + (h >> 2); 
            }
            return h;
        }
    };
    

    std::ifstream in(filename);
    std::string token;

	while(token != "POINTS")
    {
		in >> token;
        if (in.eof())
            break;
    }
    
    int n_points;
	in >> n_points;

	in >> token; //double

	V.resize(n_points, 3);

	for(int i=0; i<n_points; i++)
        for (int  j = 0; j < 3; j++)
            in >> V(i, j);
        
    while(token != "CELLS")
    {
		in >> token;
        if (in.eof())
            break;
    }
    int n_cells, n_entries;
	in >> n_cells;
	in >> n_entries;
    
    int cell_type;
	std::vector<IV4> tets;
	for(int i=0; i<n_cells; ++i)
	{
		in >> cell_type;

		if(cell_type == 4)
		{
            IV4 tet;
			for(int j = 0; j < 4; j++)
				in >> tet[j];
            tets.push_back(tet);
		}
		else
		{
			// type 1 2
			for(int j = 0; j < cell_type; j++)
				in >> token;
		}
	}
    int n_tets = tets.size();
    T.resize(n_tets, 4);
    for (int i = 0; i < n_tets; i++)
        T.row(i) = tets[i];
    
    in.close();

    std::unordered_set<IV3, IV3Hash> surface;
    for (IV4& tet : tets)
    {
        for (IV3 face : { IV3(tet[0], tet[1], tet[2]), 
                    IV3(tet[0], tet[2], tet[3]), 
                    IV3(tet[0], tet[3], tet[1]),
                    IV3(tet[1], tet[3], tet[2])})
        {
            int previous_size = surface.size();
            IV3 fi = IV3(face[0], face[1], face[2]);
            bool find = surface.find(fi) != surface.end();
            if (find) surface.erase(fi);
            fi = IV3(face[0], face[2], face[1]);
            find = surface.find(fi) != surface.end();
            if (find) surface.erase(fi);
            fi = IV3(face[1], face[0], face[2]);
            find = surface.find(fi) != surface.end();
            if (find) surface.erase(fi);
            fi = IV3(face[1], face[2], face[0]);
            find = surface.find(fi) != surface.end();
            if (find) surface.erase(fi);
            fi = IV3(face[2], face[0], face[1]);
            find = surface.find(fi) != surface.end();
            if (find) surface.erase(fi);
            fi = IV3(face[2], face[1], face[0]);
            find = surface.find(fi) != surface.end();
            if (find) surface.erase(fi);

            if (previous_size == surface.size())
            {
                surface.insert(IV3(face[2], face[1], face[0]));
            }
            
        }
    }
    
    F.resize(surface.size(), 3);
    int face_cnt = 0;
    for (const auto& face : surface)
    {
        F.row(face_cnt++) = face;
    }
}

void loadMeshFromVTKFile(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    using TV3 = Vector<T, 3>;
    using IV3 = Vector<int, 3>;

    

    std::ifstream in(filename);
    std::string token;

	while(token != "POINTS")
    {
		in >> token;
        if (in.eof())
            break;
    }
    
    int n_points;
	in >> n_points;

	in >> token; //double

	V.resize(n_points, 3);

	for(int i=0; i<n_points; i++)
        for (int  j = 0; j < 3; j++)
            in >> V(i, j);
        
    while(token != "CELLS")
    {
		in >> token;
        if (in.eof())
            break;
    }
    int n_cells, n_entries;
	in >> n_cells;
	in >> n_entries;
    
    int cell_type;
	std::vector<Vector<int, 3>> faces;
	for(int i=0; i<n_cells; ++i)
	{
		in >> cell_type;

		if(cell_type == 3)
		{
            IV3 face;
			for(int j = 0; j < 3; j++)
				in >> face[j];
            faces.push_back(face);
		}
		else
		{
			// type 1 2
			for(int j = 0; j < cell_type; j++)
				in >> token;
		}
	}
    int n_faces = faces.size();
    F.resize(n_faces, 3);
    for (int i = 0; i < n_faces; i++)
    {
        TV3 ei = V.row(faces[i][0]) - V.row(faces[i][1]);
        TV3 ej = V.row(faces[i][2]) - V.row(faces[i][1]);
        if (ej.cross(ei).dot(TV3(0, 0, 1)) < 0)
            F.row(i) = IV3(faces[i][1], faces[i][0], faces[i][2]);
        else
            F.row(i) = faces[i];
    }
    in.close();
}