#ifndef DIJETCORE_LIB_RANDOM_PT_DISTRIBUTION_H
#define DIJETCORE_LIB_RANDOM_PT_DISTRIBUTION_H


// an implementation of a PDF of x*exp(-x/T) that can be used with
// the RNGs in <random>, for a simple background generator. It uses
// the simpson rule for integration - fast, but requires more bins for
// accuracy. Function is approximated by a parabola at each point, when
// generating the random number.

#include <limits>
#include <ostream>
#include <istream>
#include <random>
#include <vector>
#include <math.h>
#include <algorithm>
#include <functional>

namespace dijetcore {
  
  template <class _RealType = double>
  class pt_distribution {
    
    typedef _RealType result_type;
    
    class param_type {
      result_type __T_;
      size_t __N_;
      result_type __Min_;
      result_type __Max_;
      std::vector<result_type> __I_;
      std::vector<result_type> __A_;
      std::vector<result_type> __B_;
      std::vector<result_type> __C_;
      bool __U_;
    public:
      typedef pt_distribution distribution_type;
      
      inline explicit param_type(result_type  _t, size_t __N, result_type __Min,
                                 result_type __Max)
      : __T_(_t), __N_(__N), __Min_(__Min), __Max_(__Max), __I_(__N_+1),
      __A_(__N_+1), __B_(__N_), __C_(__N_), __U_(false) {}
      
      inline result_type temperature() const {return __T_;}
      inline size_t bins() const {return __N_;}
      inline result_type min() const {return __Min_;}
      inline result_type max() const {return __Max_;}
      void integrate(const std::function<double(double)>& f)  {
        // Use Simpson's rule because from testing it is "good enough"
        if (!__U_) {
          result_type *x = new result_type[__N_ + 1];
          x[0] = __Min_;
          result_type dx = (__Max_ - __Min_) / __N_;
          __I_[0] = 0;
          __A_[__N_] = 0;
          for (size_t i = 0; i < __N_; ++i) {
            x[i+1] = x[i] + dx;
            __I_[i+1] = __I_[i] + ((x[i+1] - x[i])/6.0) *
            (f(x[i]) + 4 * f((x[i] + x[i+1])/2.0) +
             f(x[i+1]));
            
          }
          result_type total = __I_[__N_];
          for (int i = 1; i <= __N_; ++i) {
            __I_[i] /= __I_[__N_];
          }
          result_type x0, r1, r2, r3;
          for (size_t i = 0; i < __N_; ++i) {
            x0 = x[i];
            r1 = (dx * 0.5)/6.0 * (f(x[i]) + 4 * f((x[i] + x[i+1])/2)
                                   + f(x[i+1]))/total;
            r2 = __I_[i+1] - __I_[i];
            r3 = 2*r2 - 4*r1;
            if (fabs(r3) > 1e-8) __C_[i] = r3 / (dx * dx);
            else __C_[i] = 0;
            __A_[i] = x0;
            __B_[i] = r2 / dx - __C_[i] * dx;
            __C_[i] *= 2.0;
          }
          
          delete [] x;
          __U_ = true;
        }
      }
      
      friend inline bool operator==(const param_type& _x, const param_type& _y) {
        return _x.__T_ == _y.__T_ && _x.__N_ == _y.__N_ && _x.__Min_ == _y.__Min
        && _x.__Max_ && _y.__Max && _x.__I_ == _y.__I_ && _x.__A_ == _y.__A_
        && _x.__B_ == _y.__B_ && _x.__C_ == _y.__C_ && _x.__U_ == _y.__U_ &&
        (_x.__U_ || _x.__U_ == _y.__U_);
      }
      friend inline bool operator!=(const param_type& _x, const param_type& _y) {
        return !(_x == _y);
      }
      
      friend pt_distribution<result_type>;
    };
    
  private:
    param_type __p_;
    
  public:
    
    // constructors
    inline explicit pt_distribution(result_type __temp = 0.291, size_t __N = 500,
                                    result_type __Min = 0, result_type __Max = 5)
    : __p_(__temp, __N, __Min, __Max) {};
    inline explicit pt_distribution(param_type _p)
    : __p_(_p) {};
    
    // generator
    template<class _RNG>
    inline result_type operator()(_RNG& __g) {
      return (*this)(__g, __p_);
    }
    template <class _RNG>
    result_type operator()(_RNG& __g, param_type& __p);
    
    // numeric properties of the distribution
    inline result_type temperature() const {return __p_.__temperature__;}
    
    inline param_type param() const {return __p_;}
    inline void param(const param_type& __p) {__p_ = __p;}
    
    inline result_type eval(const result_type& x) {return x * exp(-x/__p_.__T_);}
    
    result_type min() const {return __p_.min();}
    result_type max() const {return __p_.max();}
    
    friend inline bool operator==(const pt_distribution& _x,
                                  const pt_distribution& _y) {
      return _x.__p_ == _y.__p_;
    }
    friend inline bool operator!=(const pt_distribution& _x,
                                  const pt_distribution& _y) {
      return !(_x == _y);
    }
  };
  
  template<class _RealType>
  template<class _RNG>
  _RealType pt_distribution<_RealType>::operator()(_RNG& __g, param_type& p) {
    
    if (!p.__U_) {
      auto fp = std::bind(&pt_distribution<result_type>::eval, this,
                          std::placeholders::_1);
      p.integrate(fp);
    }
    
    // get random number
    std::uniform_real_distribution<_RealType> _Uniform(0, 1);
    _RealType __Rnd = _Uniform(__g);
    
    size_t _B = std::distance(p.__I_.begin(),
                              std::lower_bound(p.__I_.begin(), p.__I_.end(), __Rnd));
    if (_B > p.__N_)
      _B = p.__N_;
    __Rnd -= p.__I_[_B];
    
    _RealType y;
    if (p.__C_[_B] != 0) {
      y = (-p.__B_[_B] + sqrt(p.__B_[_B] * p.__B_[_B] + 2 * p.__C_[_B] * __Rnd)) / p.__C_[_B];
    }
    else {
      y = __Rnd / p.__B_[_B];
    }
    _RealType x = p.__A_[_B] + y;
    
    return x;
  }
  
} // namespace dijetcore

#endif // DIJETCORE_LIB_STATISTICS_PT_DISTRIBUTION_H
