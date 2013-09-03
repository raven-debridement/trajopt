#ifndef __UTILS_H__
#define __UTILS_H__

#include <math.h>
#include "collision.h"

inline double frandom() {
	return ((double) rand()) / ((double)RAND_MAX+ 1.0);
}

inline double normal() {
	double u_1 = 0;
	while (u_1 == 0) {
		u_1 = frandom();
	}
	double u_2 = 0;
	while (u_2 == 0) {
		u_2 = frandom();
	}
	return sqrt(-2*log(u_1)) * sin(2*M_PI*u_2);
}

template <size_t size>
inline Matrix<size> sampleGaussian(const Matrix<size>& mean, const Matrix<size, size>& var) {
	Matrix<size> sample;
	for (int j = 0; j < size; ++j) {
		sample[j] = normal();
	}
	Matrix<size, size> SVec, SVal;
	jacobi(var, SVec, SVal);
	for (int i = 0; i < size; ++i) {
		SVal(i,i) = sqrt(SVal(i,i));
	}
	return SVec * SVal * sample + mean;
}

inline Matrix<3,3> cpMatrix(const Matrix<3,1>& a) 
{
	Matrix<3,3> A;
	A(0,0) = 0;       A(0,1) = -a(2,0); A(0,2) = a(1,0);
	A(1,0) = a(2,0);  A(1,1) = 0;       A(1,2) = -a(0,0);
	A(2,0) = -a(1,0); A(2,1) = a(0,0);  A(2,2) = 0;

	return A;
}

inline Matrix<3,3> rotFromErr(const Matrix<3,1>& q) {
	double rr = q(0,0)*q(0,0)+q(1,0)*q(1,0)+q(2,0)*q(2,0);
	if (rr == 0) {
		return identity<3>();
	} else {
		double r = sqrt(rr);
		//printf_s("%24.24g ", cos(r) - 1);
		return cpMatrix(q * (sin(r) / r)) + identity<3>() * cos(r) + (q*~q) * ((1 - cos(r)) / rr);
		//return cpMatrix(q * (sin(r) / r)) + identity<3>() * cos(r) + (q*~q)/rr - (q*~q)*(cos(r)/rr);
	}
}

inline Matrix<3,1> errFromRot(const Matrix<3,3>& R) {
	Matrix<3,1> q;
	q(0,0) = R(2,1) - R(1,2);
	q(1,0) = R(0,2) - R(2,0);
	q(2,0) = R(1,0) - R(0,1);

	double r = sqrt(q(0,0)*q(0,0)+q(1,0)*q(1,0)+q(2,0)*q(2,0));
	double t = R(0,0) + R(1,1) + R(2,2) - 1;

	if (r == 0) {
		return zeros<3,1>();
	} else {
		return q * (atan2(r, t) / r);
	}
}

inline Matrix<4,4> transFromErr(const Matrix<6,1>& x) {
	Matrix<4,4> X;
	X = identity<4>();
	X.insert(0,0, rotFromErr(x.subMatrix<3,1>(3,0)));
	X.insert(0,3, x.subMatrix<3,1>(0,0));
	return X;
}

inline Matrix<4,1> randQuat() {
	double s = random();
	double s1 = sqrt(1 - s);
	double s2 = sqrt(s);
	double t1 = 2 * M_PI * random();
	double t2 = 2 * M_PI * random();

	Matrix<4,1> q;

	q(0,0) = sin(t1) * s1; // x
	q(1,0) = cos(t1) * s1; // y
	q(2,0) = sin(t2) * s2; // z
	q(3,0) = cos(t2) * s2; // w

	return q;
}

inline Matrix<3,3> rotFromQuat(const Matrix<4,1>& q) {
	double x = q(0,0);
	double y = q(1,0);
	double z = q(2,0);
	double w = q(3,0);
	Matrix<3,3> R;
	R(0,0) =  1 - 2*y*y - 2*z*z; R(0,1) = 2*x*y - 2*z*w;     R(0,2) = 2*x*z + 2*y*w;
	R(1,0) =  2*x*y + 2*z*w;     R(1,1) = 1 - 2*x*x - 2*z*z; R(1,2) = 2*y*z - 2*x*w;
	R(2,0) =  2*x*z - 2*y*w;     R(2,1) = 2*y*z + 2*x*w;     R(2,2) = 1 - 2*x*x - 2*y*y;
	return R;
}

inline Matrix<4,1> quatFromRot(const Matrix<3,3>& R) {
	double x = R(2,1) - R(1,2);
	double y = R(0,2) - R(2,0);
	double z = R(1,0) - R(0,1);
	double r = sqrt(x*x+y*y+z*z);
	double t = R(0,0) + R(1,1) + R(2,2);
	double angle = atan2(r,t-1);
	if (angle != 0) {
		x /= r;
		y /= r;
		z /= r;
	} else {
		x = 0;
		y = 0;
		z = 0;
	}
	Matrix<4,1> q;
	q(0,0) = sin(angle/2)*x;
	q(1,0) = sin(angle/2)*y;
	q(2,0) = sin(angle/2)*z;
	q(3,0) = cos(angle/2);

	return q;
}

struct Vertex {
  float x, y, z;
  Vertex(float vx, float vy, float vz) : x(vx),y(vy), z(vz) {};
};

// Load obj file (obstacle definitions)
inline int objLoader(char * filename, char * label, int group_id) 
{
	char line[1024];
	std::ifstream fptr (filename, std::ios::in);
	if (!fptr) {
		std::cout << std::endl;
		std::cerr << "Error: Could not open OBJ file: " << filename << std::endl;
		exit(-1);
	}

	bool success = true;

	std::vector<Vertex> verts;
	std::vector<int> vindices;

	bool has_normals = false;

	while (!fptr.eof() && success) 
	{
		fptr.getline(line, 1023);

		std::string linestr(line);
		std::string op;
		float x, y, z;

		if (linestr.empty()) {
			success = true;
			continue;
		}

		std::stringstream ss(  std::stringstream::in | std::stringstream::out);
		ss.str(line);
		ss >> op;

		// Ignore comments
		if (op.compare("#") == 0) {
			continue;
		}

		if (op.compare("v") == 0) {
			if (ss >> x >> y >> z) {
				verts.push_back(Vertex(x*0.1f, y*0.1f, z*0.1f));
				//std::cout << "Vertex: " << x << " " << y << " " << z << std::endl;
			} else {
				success = false;
			}
		}

		if (op.compare("vn") == 0) {
			if (ss >> x >> y >> z) {
				has_normals = true;
			} else {
				success = false;
			}
		}

		// TODO: Handle degenerate triangles?
		if (op.compare("f") == 0) 
		{
			// CASE: We only have vertices
			if (!has_normals) {
				int a, b, c;
				if (ss >> a >> b >> c) {
					vindices.push_back(a-1); 
					vindices.push_back(b-1); 
					vindices.push_back(c-1);
				} else { 
					success = false;
				}
			}

			// CASE: We have a full-blown .obj file
			else {
				int index[3];

				for (int i = 0; i < 3; ++i) {
					std::string group;
					ss >> group;
					// Replace slashes with whitespace
					size_t slashidx;
					slashidx = group.find_first_of('/');
					while (slashidx != std::string::npos) {
						group[slashidx] = ' ';
						slashidx = group.find_first_of('/',slashidx+1);
					}

					// Parse data
					std::stringstream gs(std::stringstream::in | std::stringstream::out);
					gs.str(group);

					int v, n;
					if (!(gs >> v >> n)) {
						success = false;
						continue;
					}

					index[i] = v - 1;
				}

				vindices.push_back(index[0]);
				vindices.push_back(index[1]);
				vindices.push_back(index[2]);

				//std::cout << "Face: " << index[0]+1 << " " << index[1]+1 << " " << index[2]+1 << std::endl;
			}
		}

		if (!success) {
			std::cout << std::endl;
			std::cerr << "Error: Bad OBJ file: " << filename << std::endl;
			exit(-1);
		}
	}

	if( fptr.is_open()) {
		fptr.close();
	}

	// Fill in vertex data
	int nv = (int)vindices.size();
	int ntris = nv/3;
	float * pts = new float[nv*3];
	int obj_id = -1;

	for(int i = 0; i < nv; ++i) 
	{
		Vertex & v = verts[vindices[i]];
		pts[3*i] = v.x;
		pts[3*i + 1] = v.y;
		pts[3*i + 2] = v.z;
	}

	verts.clear();
	vindices.clear();

  create_triangles(group_id, ntris, pts, &obj_id, label);

	//CAL_CreateTriangles(group_id, ntris, pts, &obj_id, label);

	//clean up
	delete [] pts;

	return obj_id;
}

#endif
