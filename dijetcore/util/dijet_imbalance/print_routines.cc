#include "dijetcore/util/dijet_imbalance/print_routines.h"
#include "dijetcore/lib/logging.h"

#include <iomanip>
#include <sstream>

namespace dijetcore {

std::vector<TPaveText>
GetGridLabelsForAj(unsigned nbins_x, unsigned nbins_y, double min_x,
                   double max_x, double min_y, double max_y,
                   double x_low_axis_margin, double x_high_axis_margin,
                   double y_low_axis_margin, double y_high_axis_margin,
                   double text_size_x, double text_size_y) {
  // first find physical spacing between each bin in x and y
  double useable_x_axis = 1.0 - x_low_axis_margin - x_high_axis_margin;
  double dx = useable_x_axis / (nbins_x);
  double useable_y_axis = 1.0 - y_low_axis_margin - y_high_axis_margin;
  double dy = useable_y_axis / (nbins_y);

  // define the coordinates for all points
  struct coordinates {
  //   coordinates() {
  //     x1 = 0;
  //     x2 = 0;
  //     y1 = 0;
  //     y2 = 0;
  //   }
  //   coordinates(const coordinates &c) : x1(c.x1), x2(c.x2), y1(c.y1), y2(c.y2) {}
    double x1;
    double x2;
    double y1;
    double y2;
  };
  std::vector<coordinates> x_axis_coordinates;
  std::vector<coordinates> y_axis_coordinates;

  for (int i = 0; i < nbins_x; ++i) {
    coordinates tmp;
    tmp.x1 = x_low_axis_margin + dx * (i + 1.0/2.0);
    tmp.x2 = tmp.x1 + text_size_x;
    tmp.y1 = (1.0 - y_high_axis_margin) + 0.01;
    tmp.y2 = tmp.y1 + text_size_y;
    x_axis_coordinates.push_back(tmp);
  }
  for (int i = 0; i < nbins_y; ++i) {
    coordinates tmp;
    tmp.x1 = (1.0 - x_high_axis_margin) + 0.01;
    tmp.x2 = tmp.x1 + text_size_x;
    tmp.y1 = y_low_axis_margin + dy * (i + 1.0/2.0);
    tmp.y2 = tmp.y1 + text_size_y;
    y_axis_coordinates.push_back(tmp);
  }

  // make TPaveText for each axis entry
  std::vector<TPaveText> ret;

  // x axis first
  for (int i = 0; i < nbins_x; ++i) {
    double x_bin_width = (max_x - min_x) / (nbins_x - 1);
    double val = min_x + x_bin_width * i;
    std::stringstream ss;
    ss << std::setprecision(2) << val;

    coordinates coord = x_axis_coordinates[i];
    TPaveText tmp(coord.x1, coord.y1, coord.x2, coord.y2, "NB NDC");
    tmp.SetFillStyle(0);
    tmp.SetBorderSize(0);
    TText *tmp_text = tmp.AddText(ss.str().c_str());
    tmp_text->SetTextSize(text_size_x);
    tmp_text->SetTextAngle(0);
    ret.push_back(tmp);
  }

  // y axis
  for (int i = 0; i < nbins_x; ++i) {
    double y_bin_width = (max_y - min_y) / (nbins_y - 1);
    double val = min_y + y_bin_width * i;
    std::stringstream ss;
    ss << std::setprecision(2) << val;

    coordinates coord = y_axis_coordinates[i];
    TPaveText tmp(coord.x1, coord.y1, coord.x2, coord.y2, "NB NDC");
    tmp.SetFillStyle(0);
    tmp.SetBorderSize(0);
    TText *tmp_text = tmp.AddText(ss.str().c_str());
    tmp_text->SetTextSize(text_size_x);
    tmp_text->SetTextAngle(270);
    ret.push_back(tmp);
  }
  return ret;
}

} // namespace dijetcore