#ifndef HPC_AUX_RGBCOLOR_H
#define HPC_AUX_RGBCOLOR_H 1

namespace hpc { namespace aux {

template<typename T>
struct RGBColor {
   RGBColor() : red(0), green(0), blue(0) {
   }
   RGBColor(T red, T green, T blue) :
      red(red), green(green), blue(blue) {
   }
   RGBColor<T> get_rgb() const {
      return *this;
   }
   RGBColor interpolate(const RGBColor<T>& other, double t) {
      return RGBColor(
	 (red + (other.red - red) * t),
	 (green + (other.green - green) * t),
	 (blue + (other.blue - blue) * t));
   }
   T red, green, blue;
};

} } // namespace aux, hpc

#endif
