#pragma once

#include <cstdio>
#include <vector>

#include <boost/polygon/voronoi.hpp>
#include <boost/polygon/polygon.hpp>
// using boost::polygon::voronoi_builder;
// using boost::polygon::voronoi_diagram;
// using boost::polygon::x;
// using boost::polygon::y;
// using boost::polygon::low;
// using boost::polygon::high;

using namespace boost::polygon;
using namespace cv;
using namespace std;

typedef double coordinate_type;
typedef point_data<coordinate_type> point_type;
typedef segment_data<coordinate_type> segment_type;
typedef rectangle_data<coordinate_type> rect_type;
typedef voronoi_diagram<coordinate_type> VD;
typedef VD::cell_type cell_type;
typedef VD::cell_type::source_index_type source_index_type;
typedef VD::cell_type::source_category_type source_category_type;
typedef VD::edge_type edge_type;

rect_type brect_;
vector<point_type> point_data_;
vector<segment_type> segment_data_;

struct voro_Point {
	int a;
	int b;
	voro_Point(int x, int y) : a(x), b(y) {}
};

struct Segment {
	voro_Point p0;
	voro_Point p1;
	Segment(int x1, int y1, int x2, int y2) : p0(x1, y1), p1(x2, y2) {}
};

segment_type retrieve_segment(const cell_type& cell)
{
	source_index_type index = cell.source_index() - point_data_.size();
	return segment_data_[index];
}

point_type retrieve_point(const cell_type& cell)
{
	source_index_type index = cell.source_index();
	source_category_type category = cell.source_category();
	if (category == SOURCE_CATEGORY_SINGLE_POINT)
	{
		return point_data_[index];
	}
	index -= point_data_.size();
	if (category == SOURCE_CATEGORY_SEGMENT_START_POINT)
	{
		return low(segment_data_[index]);
	}
	else
	{
		return high(segment_data_[index]);
	}
}

namespace boost {
	namespace polygon {

		template <>
		struct geometry_concept<voro_Point> {
			typedef point_concept type;
		};

		template <>
		struct point_traits<voro_Point>
		{
			typedef int coordinate_type;

			static inline coordinate_type get(const voro_Point& point, orientation_2d orient)
			{
				return (orient == HORIZONTAL) ? point.a : point.b;
			}
		};

		template <>
		struct geometry_concept<Segment> {
			typedef segment_concept type;
		};

		template <>
		struct segment_traits<Segment> {
			typedef int coordinate_type;
			typedef voro_Point point_type;

			static inline point_type get(const Segment& segment, direction_1d dir) {
				return dir.to_int() ? segment.p1 : segment.p0;
			}
		};
	}  // polygon
}  // boost
