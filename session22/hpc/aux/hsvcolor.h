#ifndef HPC_AUX_HSVCOLOR_H
#define HPC_AUX_HSVCOLOR_H 1

#include <cmath>
#include <hpc/aux/rgbcolor.h>

namespace hpc { namespace aux {

/* the conversion function is based on the code by Nan C. Schaller at
   http://www.cs.rit.edu/~ncs/color/t_convert.html */
template<typename T>
RGBColor<T> HSVtoRGB(T h, T s, T v) {
   RGBColor<T> result;
   if (s == 0) {
      // achromatic (grey)
      return RGBColor<T>(v, v, v);
   }
   h /= 60;
   int i = std::floor(h);
   T f = h - i;
   T p = v * (1 - s);
   T q = v * (1 - s * f);
   T t = v * (1 - s * (1 - f));

   switch (i) {
      case 0: return RGBColor<T>(v, t, p);
      case 1: return RGBColor<T>(q, v, p);
      case 2: return RGBColor<T>(p, v, t);
      case 3: return RGBColor<T>(p, q, v);
      case 4: return RGBColor<T>(t, p, v);
      default: return RGBColor<T>(v, p, q);
   }
}

template<typename T>
struct HSVColor {
   HSVColor() : hue(0), saturation(0), value(0) {
   }
   HSVColor(T hue, T saturation, T value) :
      hue(hue), saturation(saturation), value(value) {
   }
   RGBColor<T> get_rgb() const {
      return HSVtoRGB(hue, saturation, value);
   }
   HSVColor<T> interpolate(const HSVColor<T>& other, T t) {
      return HSVColor(
	 (other.hue == hue &&
	    other.saturation == saturation && other.value == value?
	       std::fmod(hue + (360 + other.hue - hue) * t, 360)
	    :
	       hue + (other.hue - hue) * t),
	 saturation + (other.saturation - saturation) * t,
	 value + (other.value - value) * t
      );
   }
   T hue;
   T saturation;
   T value;
};

} } // namespace aux, hpc

#endif
