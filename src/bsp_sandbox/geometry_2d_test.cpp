#include "geometry_2d.hpp"

using namespace Geometry2D;

int main() {
  {
    Point p; p << 0, 0;
    Point v; v << 0, 1;
    Point w; w << 1, 0;
    cout << "should be 0.707: " << dist_to_segment(p, v, w) << endl;
  }
  {
    Point p1; p1 << 0, 0;
    Point p2; p2 << 1, 1;
    Point p3; p3 << 1, 0;
    Point p4; p4 << 0, 1;

    Point intersection;
    cout << "should be 1: " << line_intersection(p1, p2, p3, p4, &intersection) << endl;
    cout << "should be (0.5, 0.5): " << intersection.transpose() << endl;
  }
  {
    Point p1; p1 << 1, 0;
    Point p2; p2 << 0, 1;
    Point p3; p3 << 0, 0;
    Point p4; p4 << 0, 0.5;
    cout << "should be 0: " << segment_intersects_segment(p1, p2, p3, p4) << endl;
  }
  {
    Point p1; p1 << 1, 0;
    Point p2; p2 << 0, 1;
    Point p3; p3 << 0, 0;
    Point p4; p4 << 0, 1;
    cout << "should be 1: " << segment_intersects_segment(p1, p2, p3, p4) << endl;
  }
  {
    Beam2D beam;
    beam.base = Point(0, 0);
    beam.a = Point(1, 3);
    beam.b = Point(-1, 3);
    Segment seg;
    seg.a = Point(-1, 4);
    seg.b = Point(1, 4);
    Segment seg2;
    seg2.a = Point(-0.5, 1.5);
    seg2.b = Point(0.5, 1.5);
    Segment seg3;
    seg3.a = Point(0, 4);
    seg3.b = Point(2, 2);
    Segment seg4;
    seg4.a = Point(-1, 1);
    seg4.b = Point(0, 4);
    Segment seg5;
    seg5.a = Point(-1, 0);
    seg5.b = Point(1, 0);
    cout << "initial beam: " << beam << endl;
    cout << "should be the same as initial: " << truncate_beam(beam, seg) << endl;
    cout << "should be truncated: " << truncate_beam(beam, seg2) << endl;
    cout << "should be the same as initial: " << truncate_beam(beam, seg3) << endl;
    cout << "should be truncated: " << truncate_beam(beam, seg4) << endl;
    cout << "should be shrunk to a point: " << truncate_beam(beam, seg5) << endl;
  }
  {
    Beam2D beam;
    beam.base = Point(0, 0);
    beam.a = Point(1, 3);
    beam.b = Point(-1, 3);
    Segment seg;
    seg.a = Point(-1, 1);
    seg.b = Point(0, 4);
    Segment seg2;
    seg2.a = Point(1, 1);
    seg2.b = Point(0, 4);
    vector<Segment> segs;
    segs.push_back(seg);
    segs.push_back(seg2);
    cout << "initial beam: " << beam << endl;
    cout << "should be truncated: " << truncate_beam(beam, segs) << endl;
    segs.push_back(seg);
    segs.push_back(seg2);
    cout << "should be same as above: " << truncate_beam(beam, segs) << endl;
    segs.push_back(seg);
    segs.push_back(seg2);
    cout << "should be same as above: " << truncate_beam(beam, segs) << endl;
  }
  {
    Beam2D beam;
    beam.base = Point(0, 0);
    beam.a = Point(1, 3);
    beam.b = Point(-1, 3);
    cout << "should be 1: " << inside(Point(0, 0), beam) << endl;
    cout << "should be 1: " << inside(Point(0, 1), beam) << endl;
    cout << "should be 1: " << inside(Point(0, 2), beam) << endl;
    cout << "should be 1: " << inside(Point(0, 3), beam) << endl;
    cout << "should be 0: " << inside(Point(0, 4), beam) << endl;
    cout << "should be 0: " << inside(Point(-1, 1), beam) << endl;
  }
  {
    Beam2D beam1;
    beam1.base = Point(0, 0);
    beam1.a = Point(0, 3);
    beam1.b = Point(-1, 3);
    Beam2D beam2;
    beam2.base = Point(0, 0);
    beam2.a = Point(2, 4);
    beam2.b = Point(0, 4);
    Beam2D beam3;
    beam3.base = Point(0, 0);
    beam3.a = Point(2, 3);
    beam3.b = Point(1, 3);
    vector<Beam2D> beams;
    beams.push_back(beam3);
    beams.push_back(beam2);
    beams.push_back(beam1);
    cout << "should be 1.0: " << sgndist(Point(-1, 4), beams) << endl;
    cout << "should be 1: " << inside(Point(0, 2), beams) << endl;
    cout << "should be " << dist_to_segment(Point(0, 2), beams[2].base, beams[2].b) << ": " << sgndist(Point(0, 2), beams) << endl;
  }
  {
    Beam2D beam;
    beam.base = Point(0, 0);
    beam.a = Point(1, 3);
    beam.b = Point(-1, 3);
    vector<Beam2D> beams;
    partition(beam, 2, &beams);
    cout << "should be half beam: " << beams[0] << endl;
    cout << "should be half beam: " << beams[1] << endl;
  }
  {
    cout << "cdf(standard_normal, 2.57) = " << cdf(standard_normal, 2.57) << endl;
  }
  {
    double mean, cov;
    double z = 0;
    truncate_univariate_gaussian(z, 0, 1, &mean, &cov);
    cout << "truncate standard gaussian (>=0). mean: " << mean << "; cov: " << cov << endl;
  }
  {
    Vector2d c(-1, -1);
    double d = 0;
    Vector2d cur_mean(0, 0);
    Matrix2d cur_cov = Matrix2d::Identity();
    Vector2d new_mean;
    Matrix2d new_cov;
    truncate_gaussian(c, d, cur_mean, cur_cov, &new_mean, &new_cov);
    cout << "truncated mean: " << new_mean.transpose() << endl;
    cout << "truncated covariance: \n" << new_cov << endl;
  }
  {
    Beam2D beam1;
    beam1.base = Point(-1, -1);
    beam1.a = Point(3, 0);
    beam1.b = Point(0, 0);
    Beam2D beam2;
    beam2.base = Point(-1, -1);
    beam2.a = Point(0, 0);
    beam2.b = Point(0, 3);
    vector<Beam2D> beams;
    beams.push_back(beam1);
    beams.push_back(beam2);
    Vector2d cur_mean(0, 0);
    Matrix2d cur_cov = Matrix2d::Identity();
    Vector2d new_mean;
    Matrix2d new_cov;
    truncate_belief(beams, cur_mean, cur_cov, &new_mean, &new_cov);
    cout << "truncated mean: " << new_mean.transpose() << endl;
    cout << "truncated covariance: \n" << new_cov << endl;
  }
  return 0;
}
