#ifndef HPC_MATVEC_MATRIX2PIXBUF_H
#define HPC_MATVEC_MATRIX2PIXBUF_H

#include <gdk-pixbuf/gdk-pixbuf.h>
#include <hpc/aux/rgbcolor.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/isgematrix.h>

namespace hpc { namespace matvec {

/*
   generate colorized pixmap from A:

      y
      ^
      |   A(0, A.numCols-1)    ...      A(A.numRows-1, A.numCols-1)
      |          .                                   .
      |          .                                   .
      |          .                                   .
      |        A(0, 0)         ...           A(A.numRows-1, 0)
      +---------------------------------------------------------------> x

   whereas the first index from 0 to A.numRows-1 is associated with the x axis
   and whereas the second index from 0 to A.numCols-1 is tied to the y axis

   the scale factor creates larger blocks of scale x scale pixels
*/
template <typename MA, typename Colorizer>
typename std::enable_if<IsGeMatrix<MA>::value,
         GdkPixbuf*>::type
create_pixbuf(MA& A, Colorizer colorizer, unsigned int scale = 1) {
   int width = A.numRows * scale;
   int height = A.numCols * scale;
   GdkPixbuf* pixbuf = gdk_pixbuf_new(GDK_COLORSPACE_RGB, FALSE,
      8, width, height);
   int n_channels = gdk_pixbuf_get_n_channels(pixbuf);
   int rowstride = gdk_pixbuf_get_rowstride(pixbuf);
   guchar* pixels = gdk_pixbuf_get_pixels(pixbuf);
   if (!pixels) return nullptr;

   for (std::size_t i = 0; i < A.numCols; ++i) {
      for (std::size_t j = 0; j < A.numRows; ++j) {
	 auto color = colorizer(A(i, j));
	 auto rgbcolor = color.get_rgb();
	 for (std::size_t scale_i = 0; scale_i < scale; ++scale_i) {
	    for (std::size_t scale_j = 0; scale_j < scale; ++scale_j) {
	       int i_ = height - 1 - (j * scale + scale_j);
	       int j_ = i * scale + scale_i;
	       assert(i_ >= 0 && i_ < height);
	       assert(j_ >= 0 && j_ < width);
	       guchar* pp = pixels + i_ * rowstride + j_ * n_channels;
	       pp[0] = rgbcolor.red * 255;
	       pp[1] = rgbcolor.green * 255;
	       pp[2] = rgbcolor.blue * 255;
	    }
	 }
      }
   }
   return pixbuf;
}

} } // namespaces matvec, hpc

#endif
