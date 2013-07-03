#pragma once
#include "bsp/bsp.hpp"

namespace Geometry2D {
  const double eps = 1e-5;
  typedef Vector2d Point;

  // Assume that the points {base, a, b} are ordered counterclockwise
  struct Beam2D {
    Point base;
    Point a;
    Point b;
    friend ostream& operator<< (ostream &out, const Beam2D& beam);
  };

  struct Segment {
    Point a;
    Point b;
  };

  double dist_to_segment(const Point& p, const Point& v, const Point& w);

  bool line_intersection(const Point& p1, const Point& p2, const Point& p3, const Point& p4, Point* intersection);

  bool segment_intersects_segment(const Point& p1, const Point& p2, const Point& p3, const Point& p4);

  Beam2D truncate_beam(const Beam2D& beam, const Segment& segment);

  Beam2D truncate_beam(const Beam2D& beam, const vector<Segment>& segments);

  void truncate_beams(const vector<Beam2D> beams, const vector<Segment>& segments, vector<Beam2D>* truncated_beams);

  bool inside(const Point& p, const Beam2D& beam);

  bool inside(const Point& p, const vector<Beam2D>& beams);

  double sgndist(const Point& p, const Beam2D& beam);

  double sgndist(const Point& p, const vector<Beam2D>& beams);

  void partition(const Beam2D& beam, int nparts, vector<Beam2D>* beams);
    
}
