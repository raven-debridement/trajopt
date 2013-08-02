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

  template<class T1, class T2>
  struct is_same_scalar_type {
    typedef typename MatrixTraits<T1>::scalar_type scalar_type_1;
    typedef typename MatrixTraits<T2>::scalar_type scalar_type_2;
    const static bool value = std::is_same<scalar_type_1, scalar_type_2>::value;
  };

  template<class T1>
  struct is_vector {
    const static bool value = MatrixTraits<T1>::cols == 1;
  };

}
