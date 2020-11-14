#include <iostream>
#include <vector>
#include <algorithm>
#include <climits>
#include <sstream>
#include <stdint.h>
#include <queue>
#include <functional>
#include <array>
#include <memory>

using std::cout;
using std::cin;
using std::endl;

#define llint int64_t // uint64_t

#define DEBUG 1

#if (DEBUG == 1)
#define disp(s, x) std::cout << s << " " << x << std::endl;
#define msg(s) std::cout << s << std::endl;
#else
#define disp(s, x) ;
#define msg(s) ;
#endif

#define print(x) std::cout << x << "\n";

// inplace map for (int64_t array -> int64_t array)
void map(llint *arr, int N, llint(*func)(llint&)) {
    for (int i=0; i < N; i++) {
        arr[i] = func(arr[i]);
    }
}

// inplace filter for (int64_t -> bool)
// returns N: new size of array
int filter(llint *arr, int N, bool(*func)(llint&)) {
    int j = 0;
    int newN = N;
    for (int i=0; i<N; i++) {
        if (func(arr[i])) {
            arr[j] = arr[i];
            j++;
        } else {
            newN--;
        }
    }
    return newN;
}

// accumulates value over array with (int64_t, int64_t -> int64_t)
llint reduce(llint *arr, int N, llint(*func)(llint&, llint&)) {
    llint val = func(arr[0], arr[1]);
    for (int i=2; i < N; i++) {
        val = func(val, arr[i]);
    }
    return val;
}


class Row {
public:
    explicit Row(size_t r) : row(r) {};
    size_t operator()() const {return row;}
private:
    size_t row;
};

class Col {
public:
    explicit Col( size_t c) : col(c) {};
    size_t operator()() const { return col; }
private:
    size_t col;
};

template <typename T>
class Matrix {
    // reference: http://www.cplusplus.com/forum/articles/17108/
public:
    Matrix(Row r, Col c) : rows_(r()), cols_(c()), arr_(0) {
        if (cols_ > 0 && rows_ > 0) {
            ref_ = new int;
            *ref_ = 0;
            arr_ = new T[rows_*cols_];
        }
    }

    Matrix(Col c, Row r) : cols_(c()), rows_(r()), arr_(0) {
        if (cols_ > 0 && rows_ > 0) {
            ref_ = new int;
            *ref_ = 0;
            arr_ = new T[rows_*cols_];
        }
    }

    Matrix(unsigned row, unsigned col) : rows_(row), cols_(col), arr_(0) {
         if (cols_ > 0 && rows_ > 0) {
            ref_ = new int;
            *ref_ = 0;
            arr_ = new T[rows_*cols_];
         }
     }

    ~Matrix() {
        if (*ref_ <= 0) {
            delete[] arr_;
            delete ref_;
        } else {
            (*ref_)--;
        }
    }

    void Fill(T val) {
        for (size_t i = 0; i < rows_; i++) {
            for (size_t j = 0; j < cols_; j++) {
                arr_[i*cols_ + j] = val;
            }
        }
    }

    void Print() {
        for (size_t i=0; i<rows_; i++) {
            for (size_t j=0; j<cols_; j++) {
                std::cout << arr_[i*cols_+j] << " ";
            }
            std::cout << "\n";
        }
    }

    size_t cols() const { return cols_; }
    size_t rows() const { return rows_; }

    // matrix access with parenthesis, i.e. M(row, col)
    const T& operator() (unsigned x, unsigned y) const { return arr_[y*cols_ + x]; }
    T& operator() (unsigned row, unsigned col) { return arr_[row*cols_ + col]; }
    T& operator() (Col x, Row y) { return arr_[y()*cols_ + x()]; }
    T& operator() (Row y, Col x) { return arr_[y()*cols_ + x()]; }


    Matrix<T> clone() {
        Matrix<T> B(rows_, cols_);
        for (size_t i = 0; i<rows_; i++) {
            for (size_t j=0; j<cols_; j++) {
                B.arr_[i*cols_ + j] = arr_[i*cols_ + j];
            }
        }
        // (*ref_)++;
        return B;
    }

    void copy(Matrix<T>& B) {
        if (B.rows() == rows_ && B.cols() == cols_) {
            for (size_t i = 0; i<rows_; i++) {
                for (size_t j=0; j<cols_; j++) {
                    arr_[i*cols_ + j] = B.arr_[i*cols_ + j];
                }
            }
        }
    }

    // in place matrix addition with broadcasting
    void add(Matrix<T>& B) {
        if ((B.cols() == cols_) && (B.rows() == rows_)) {
            for (size_t i = 0; i<rows_; i++) {
                for (size_t j=0; j<cols_; j++) {
                    arr_[i*cols_ + j] += B.arr_[i*cols_ + j];
                }
            }
        } else if (B.rows() == 1) {
            for (size_t i = 0; i<rows_; i++) {
                for (size_t j=0; j<cols_; j++) {
                    arr_[i*cols_ + j] += B.arr_[i];
                }
            }
        } else if (B.cols() == 1) {
            for (size_t i = 0; i<rows_; i++) {
                for (size_t j=0; j<cols_; j++) {
                    arr_[i*cols_ + j] += B.arr_[j];
                }
            }
        }
    }
    
    // in place matrix subtraction with broadcasting
    void sub(Matrix<T>& B) {
        if ((B.cols() == cols_) && (B.rows() == rows_)) {
            for (size_t i = 0; i<rows_; i++) {
                for (size_t j=0; j<cols_; j++) {
                    arr_[i*cols_ + j] -= B.arr_[i*cols_ + j];
                }
            }
        } else if (B.rows() == 1) {
            for (size_t i = 0; i<rows_; i++) {
                for (size_t j=0; j<cols_; j++) {
                    arr_[i*cols_ + j] -= B.arr_[i];
                }
            }
        } else if (B.cols() == 1) {
            for (size_t i = 0; i<rows_; i++) {
                for (size_t j=0; j<cols_; j++) {
                    arr_[i*cols_ + j] -= B.arr_[j];
                }
            }
        }
    }

    // in place addition with a number
    void add(T& num) {
        for (size_t i=0; i< rows_; i++) {
            for (size_t j=0; j<cols_; j++) {
                arr_[i*cols_ + j] += num;
            }
        }
    }

    // in place multiplication by a number
    void mult(T num) {
        for (size_t i=0; i< rows_; i++) {
            for (size_t j=0; j<cols_; j++) {
                arr_[i*cols_ + j] *= num;
            }
        }
    }

    // multiplication by a number (not in place)
    Matrix<T> operator* (T x) {
        // declares array tmp[rows_*cols_];
        Matrix<T> B(rows_, cols_);

        for (size_t i=0; i< rows_; i++) {
            for (size_t j=0; j<cols_; j++) {
                B.arr_[i*cols_ + j] = arr_[i*cols_ + j] * x;
            }
        }
        // prevent a double free of B.arr_
        (*B.ref_) += 1;
        return B;
    }

    // matrix multiplication
    Matrix<T> operator* (Matrix<T>& B) { return this->dot(B); }

    // matrix multiplication (alternate syntax to A*B)
    Matrix<T> dot(Matrix<T>& B) {
        // TODO: implement strassen multiplication
        // C11 = M1 + M4 - M5 + M7
        // C12 = M3 + M5
        // C21 = M2 + M4
        // C22 = M1 - M2 + M3 + M6

        // M1 = (A11 + A22)*(B11 + B22)
        // M2 = (A21 + A22)*B11
        // M3 = A11*(B12 - B22)
        // M4 = A22*(B21 - B11)
        // M5 = (A11 + A12)*B22
        // M6 = (A21 - A11)*(B11 + B12)
        // M7 = (A12 - A22)*(B21 + B22)

        // Matrix<T> C11(Row(rows_), Col(cols_));
        // Matrix<T> C12(Row(rows_), Col(cols_));
        // Matrix<T> C21(Row(rows_), Col(cols_));
        // Matrix<T> C22(Row(rows_), Col(cols_));
        // Matrix<T> M1(Row(rows_), Col(cols_));
        // Matrix<T> M2(Row(rows_), Col(cols_));
        // Matrix<T> M3(Row(rows_), Col(cols_));
        // Matrix<T> M4(Row(rows_), Col(cols_));
        // Matrix<T> M5(Row(rows_), Col(cols_));
        // Matrix<T> M6(Row(rows_), Col(cols_));
        // Matrix<T> M7(Row(rows_), Col(cols_));

        Matrix<T> C(Row(rows_), Col(B.cols()));
        T x;
        for (size_t r1=0; r1<rows_; r1++) {
            for (size_t c2=0; c2<B.cols(); c2++) {
                x = 0;
                for (size_t c1=0; c1<cols_; c1++) {
                    x += arr_[r1*cols_ + c1] * B.arr_[c1*B.cols() + c2];
                }
                C(Row(r1), Col(c2)) = x;
            }
        }
        return C;
    }

    // Increment reference
    // prevents destructor from freeing pointers 1 time
    // use ONLY when you need to use the matrix outside of the scope it was declared
    void incRef() { (*ref_)++; }

    // void operator= (Matrix<T>& B) {
    //     this->copy(B);
    // }

private:
    size_t cols_;
    size_t rows_;
    T* arr_;
    int* ref_;

    // prevents unwanted copying
    // Matrix(const Matrix<T>&);
    // Matrix& operator = (const Matrix<T>&);
};

template <typename T, typename T2>
Matrix<T> operator* (T2 num, Matrix<T>& A) { return A*num; }

template <typename T, typename T2>
Matrix<T> operator- (T2 num, Matrix<T>& A) {
    T val = (T) num;
    size_t rows = A.rows();
    size_t cols = A.cols();
    Matrix<T> B(rows, cols);
    B.copy(A);
    B.mult(-1);
    B.add(val);
    B.incRef();
    return B;
}

template <typename T, typename T2>
Matrix<T> operator- (Matrix<T>& A, T2 num) {
    T val = (T) -1*num;
    size_t rows = A.rows();
    size_t cols = A.cols();
    Matrix<T> B(rows, cols);
    B.copy(A);
    B.add(val);
    B.incRef();
    return B;
}

template <typename T, typename T2>
Matrix<T> operator+ (Matrix<T>& A, T2 num) {
    T val = (T) num;
    size_t rows = A.rows();
    size_t cols = A.cols();
    Matrix<T> B(rows, cols);
    B.copy(A);
    B.add(val);
    B.incRef();
    return B;
}

template <typename T, typename T2>
Matrix<T> operator+ (T2 num, Matrix<T>& A) {
    T val = (T) num;
    size_t rows = A.rows();
    size_t cols = A.cols();
    Matrix<T> B(rows, cols);
    B.copy(A);
    B.add(val);
    B.incRef();
    return B;
}

template <typename T, typename T2>
Matrix<double> operator/ (Matrix<T2>& A, T num) {
    size_t rows = A.rows();
    size_t cols = A.cols();
    Matrix<double> B(rows, cols);
    B.copy(A);
    for (size_t i=0; i<rows; i++) {
        for(size_t j=0; j<cols; j++) {
            B(Row(rows), Col(cols)) /= num;
        }
    }
    B.incRef();
    return B;
}


void printarr(llint *arr, int N) {
    for (int i=0; i<N; i++) {
        cout << arr[i] << "\n";
    }
}

// Copies A to B
void cpy(llint* A, llint* B, int N) {
    for (int i = 0; i < N; i++) {
        B[i] = A[i];
    }
}

// int main(int argc, char **argv) {
int main() {
    llint N = 4;
    llint A[4] = {1,2,3,4};
    print("A");
    printarr(A, N);

    // some sample functions
    auto times2 = [](auto& x) { return x*2; };
    auto add = [](auto& x, auto& y) {return x+y;};
    auto isOdd = [](auto& x) -> bool {return x%2 == 1; };
    // auto isEven = [](auto& x) -> bool {return x % 2 == 0;};

    print("\nfilter: odds")
    N = filter(A, N, isOdd);
    printarr(A, N);

    print("\nmap: times 2")
    map(A, N, times2);

    printarr(A, N);
    print("\nreduce: with sum")

    llint sum = reduce(A, N, add);
    print(sum);

    Matrix<llint> V1(Row(2), Col(3));
    Matrix<llint> V2(Row(3), Col(2));

    V1.Fill(5);
    V2.Fill(3);

    print("\nV1");
    V1.Print();
    print("\nV2");
    V2.Print();

    print("\nV1 * V2");

    Matrix<llint> V3 = V1*V2;
    V3.Print();

    print("\nV1.dot(V2)");
    Matrix<llint> V4 = V1.dot(V2);
    V4.Print();

    print("\nMatrix accessing");
    V3(1,1) = -2;
    V3(Col(1), Row(0)) = 8;
    V3(Row(1), Col(0)) *= -1;
    V3.Print();


    print("\n2*V3:\n");
    V3 = 2*V3;
    V3.Print();

    print("\n(V3) - 1");
    V3 = V3 - 1;
    V3.Print();

    print("\n5 - (V3) + 7");
    V3 = (5 - V3);
    V3 = V3 - 7;
    V3.Print();

    return 0;
}