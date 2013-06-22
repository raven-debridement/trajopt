#pragma once

namespace BSP {
  template<class T> struct MatrixTraits;
  template<
    template<typename, int, int, int, int, int> class MatrixType,
    typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols
  >
  struct MatrixTraits< MatrixType<Scalar, Rows, Cols, Options, MaxRows, MaxCols> > {
    const static int rows = Rows;
    const static int cols = Cols;
    const static int options = Options;
    const static int max_rows = MaxRows;
    const static int max_cols = MaxCols;
    typedef Scalar scalar_type;
  };
}
