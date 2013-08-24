#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkOBJReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkSelectEnclosedPoints.h>
#include <string>

using namespace std;

struct PointGenerator {
  double min_v[3];
  double max_v[3];
  PointGenerator(const double min_coords[3], const double max_coords[3]) {
    memcpy(min_v, min_coords, sizeof(min_v));
    memcpy(max_v, max_coords, sizeof(max_v));
  }
  double runif(double a, double b) {
    return a + (((double) rand()) / (double) RAND_MAX) * (b - a);
  }
  void generate(double d[3]) {
    d[0] = runif(min_v[0], max_v[0]);
    d[1] = runif(min_v[1], max_v[1]);
    d[2] = runif(min_v[2], max_v[2]);
  }
};


bool point_in_mesh(double pt[], vtkSmartPointer<vtkPolyData> mesh) {
  vtkSmartPointer<vtkPoints> points = 
    vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(pt);

  vtkSmartPointer<vtkPolyData> pointsPolydata = 
    vtkSmartPointer<vtkPolyData>::New();

  pointsPolydata->SetPoints(points);
  vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPoints = 
    vtkSmartPointer<vtkSelectEnclosedPoints>::New();
  selectEnclosedPoints->SetInputData(pointsPolydata);
  selectEnclosedPoints->SetSurfaceData(mesh);
  selectEnclosedPoints->Update();

  return selectEnclosedPoints->IsInside(0);
}

vtkSmartPointer<vtkPolyData> get_mesh_by_filename(string filename) {
  vtkSmartPointer<vtkOBJReader> reader =
    vtkSmartPointer<vtkOBJReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  vtkSmartPointer<vtkPolyData> mesh = reader->GetOutput();

  mesh->GetPolys()->InitTraversal();

  return mesh;
}

int main(int argc, char *argv[]) {
  // Parse command line arguments
  //
  string included_filename = "./objs/prostate.obj";
  vector<string> excluded_filenames;
  excluded_filenames.push_back("./objs/arteriesBaseLeft.obj");
  excluded_filenames.push_back("./objs/arteriesTubeLeft.obj");
  excluded_filenames.push_back("./objs/bladder.obj");
  excluded_filenames.push_back("./objs/cowpersGlands1.obj");
  excluded_filenames.push_back("./objs/dermis.obj");
  excluded_filenames.push_back("./objs/glan.obj");
  excluded_filenames.push_back("./objs/pelvis.obj");
  excluded_filenames.push_back("./objs/urethra.obj");
  excluded_filenames.push_back("./objs/arteriesBaseRight.obj");
  excluded_filenames.push_back("./objs/arteriesTubeRight.obj");
  excluded_filenames.push_back("./objs/corpusCavernosum.obj");
  excluded_filenames.push_back("./objs/cowpersGlands2.obj");
  excluded_filenames.push_back("./objs/epidermis.obj");
  excluded_filenames.push_back("./objs/hypodermis.obj");

  if(argc != 2) {
    cout << "Usage: " << argv[0] << " PointCount" << endl;
    return EXIT_FAILURE;
  }
 
  vtkSmartPointer<vtkPolyData> included_mesh = get_mesh_by_filename(included_filename);


  vector< vtkSmartPointer<vtkPolyData> > excluded_meshes;
  for (int i = 0; i < excluded_filenames.size(); ++i) {
    excluded_meshes.push_back(get_mesh_by_filename(excluded_filenames[i]));
  }

  vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkCell> cell = included_mesh->GetCell(1);
  double min_coords[3], max_coords[3];
  double read_pt[3];

  included_mesh->GetCell(0)->GetPoints()->GetPoint(0, min_coords);
  included_mesh->GetCell(0)->GetPoints()->GetPoint(0, max_coords);

  
  for (int i = 0; i < included_mesh->GetNumberOfCells(); ++i) {
    vtkSmartPointer<vtkPoints> points = included_mesh->GetCell(i)->GetPoints();
    for (int j = 0; j < points->GetNumberOfPoints(); ++j) {
      points->GetPoint(j, read_pt);
      for (int k = 0; k < 3; ++k) {
        if (read_pt[k] < min_coords[k]) {
          min_coords[k] = read_pt[k];
        }
        if (read_pt[k] > max_coords[k]) {
          max_coords[k] = read_pt[k];
        }
      }
    }
  }

  PointGenerator g(min_coords, max_coords);

  int exp_cnt = atoi(argv[1]);
  int npts = exp_cnt * 20;
  double **pt = new double*[npts];

  for (int i = 0; i < npts; ++i) {
    pt[i] = new double[3];
    g.generate(pt[i]);
  }

  int cnt = 0;
  for (int i = 0; i < npts; ++i) {
    //double enclosed_pts[28][3];
    //for (int k = 0; k <= 27; ++k)
    //  for (int j = 0; j < 3; ++j) enclosed_pts[k][j] = pt[i][j];
    //int k = 1;
    //for (int d1 = -1; d1 <= 1; ++d1) {
    //  for (int d2 = -1; d2 <= 1; ++d2) {
    //    for (int d3 = -1; d3 <= 1; ++d3) {
    //      enclosed_pts[k][0] += 0.25 * d1;
    //      enclosed_pts[k][1] += 0.25 * d2;
    //      enclosed_pts[k][2] += 0.25 * d3;
    //      ++k;
    //    }
    //  }
    //}
    //bool accepted = true;
    //for (int j = 0; j <= 27; ++j) {
    //  if (!point_in_mesh(enclosed_pts[j], included_mesh)) {
    //    accepted = false;
    //    break;
    //  }
    //}
    //
    //for (int k = 0; k <= 27; ++k) {
    //  for (int j = 0; j < excluded_meshes.size(); ++j) {
    //    if (point_in_mesh(enclosed_pts[k], excluded_meshes[j])) {
    //      accepted = false;
    //      break;
    //    }
    //  }
    //  if (!accepted) break;
    //}
    if (point_in_mesh(pt[i], included_mesh)) {
      cout << pt[i][0] << " " << pt[i][1] << " " << pt[i][2] << endl;
      ++cnt;
      if (cnt == exp_cnt) {
        break;
      }
    }
  }

  for (int i = 0; i < npts; ++i) {
    delete [] pt[i];
  }
  delete [] pt;
 
  return EXIT_SUCCESS;
}
