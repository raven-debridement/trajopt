#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace std;

struct Vertex {
  double x, y, z;
  Vertex(): x(0.0), y(0.0), z(0.0) {}
  Vertex(double x, double y, double z): x(x), y(y), z(z) {}
  Vertex& operator=(const Vertex& other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    return *this;
  }
  Vertex operator-() const {
    return Vertex(-this->x, -this->y, -this->z);
  }
  friend Vertex operator+(const Vertex& l, const Vertex& r) {
    return Vertex(l.x + r.x, l.y + r.y, l.z + r.z);
  }
  friend Vertex operator-(const Vertex& l, const Vertex& r) {
    return Vertex(l.x - r.x, l.y - r.y, l.z - r.z);
  }
  friend double operator*(const Vertex& l, const Vertex& r) {
    return l.x*r.x + l.y*r.y + l.z*r.z;
  }
  friend Vertex operator*(double c, const Vertex& v) {
    return Vertex(c*v.x, c*v.y, c*v.z);
  }
  friend Vertex operator*(const Vertex& v, double c) {
    return Vertex(c*v.x, c*v.y, c*v.z);
  }
  friend Vertex operator^(const Vertex& l, const Vertex& r) {
    return Vertex(l.y*r.z - r.y*l.z, - l.x*r.z + r.x*l.z, l.x*r.y - l.y*l.x);
  }
  friend ostream& operator<<(ostream& out, const Vertex& v) {
    out << "v " << v.x << " " << v.y << " " << v.z;
    return out;
  }
};

struct Face {
  int a, b, c;
  Face(): a(0), b(0), c(0) {}
  Face(int a, int b, int c): a(a), b(b), c(c) {}
  friend ostream& operator<<(ostream& out, const Face& f) {
    out << "f " << f.a+1 << " " << f.b+1 << " " << f.c+1;
    return out;
  }
};

int main(int argc, char**argv) {
  ifstream objfile;
  ofstream outfile;
  vector<Vertex> vertices;
  vector<Face> faces;
  double x, y, z;
  int a, b, c;
  if (argc != 4) {
    cout << "Usage: " << argv[0] << " InFilename(.obj) OutFilename(.obj) NewScale" << endl;
    return 1;
  }
  double scale = atof(argv[3]);
  objfile.open(argv[1], ios::in);
  while (!objfile.eof()) {
    string line;
    getline(objfile, line, '\r');
    stringstream ss;
    char head;
    ss << line;
    ss >> head;
    if (head == 'v') {
      ss >> x >> y >> z;
      vertices.push_back(Vertex(x, y, z));
    } else if (head == 'f') {
      ss >> a >> b >> c;
      faces.push_back(Face(a-1, b-1, c-1));
    }
  }
  objfile.close();
  outfile.open(argv[2], ios::out);
  for (vector<Vertex>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
    outfile << scale*(*it) << endl;
  }
  for (vector<Face>::iterator it = faces.begin(); it != faces.end(); ++it) {
    outfile << *it << endl;
  }
  return 0;
}
