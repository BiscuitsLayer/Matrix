#pragma once

//  SYSTEM
#include <algorithm>

template <typename T>
class MatrixBuffer {
    private:
        /* empty */
    protected:
        //  DATA
        int size_ = 0, used_ = 0;
        T* data_ = nullptr;

        //  CTOR
        MatrixBuffer (int size = 0):
            size_ (size),
            used_ (0),
            data_ (nullptr)
            {
                int memorySize = size * sizeof (T);
                data_ = (size == 0 ? nullptr : static_cast <T*> (::operator new [] (memorySize)));
            }

        //  METHOD
        void Swap (MatrixBuffer& rhs) {
            std::swap (this->size_, rhs.size_);
            std::swap (this->used_, rhs.used_);
            std::swap (this->data_, rhs.data_);
        }
        
        //  DTOR
        ~MatrixBuffer () {
            for (int i = 0; i < used_; ++i) {
                data_[i].~T ();
            }
            operator delete [] (data_);
        }
    public:
        //  CTOR
        MatrixBuffer (const MatrixBuffer& rhs) = delete;

        //  OVERLOADED OPERATOR
        MatrixBuffer& operator = (const MatrixBuffer& rhs) = delete;
};