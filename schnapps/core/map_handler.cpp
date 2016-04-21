/*******************************************************************************
* SCHNApps                                                                     *
* Copyright (C) 2015, IGG Group, ICube, University of Strasbourg, France       *
*                                                                              *
* This library is free software; you can redistribute it and/or modify it      *
* under the terms of the GNU Lesser General Public License as published by the *
* Free Software Foundation; either version 2.1 of the License, or (at your     *
* option) any later version.                                                   *
*                                                                              *
* This library is distributed in the hope that it will be useful, but WITHOUT  *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  *
* for more details.                                                            *
*                                                                              *
* You should have received a copy of the GNU Lesser General Public License     *
* along with this library; if not, write to the Free Software Foundation,      *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.           *
*                                                                              *
* Web site: http://cgogn.unistra.fr/                                           *
* Contact information: cgogn@unistra.fr                                        *
*                                                                              *
*******************************************************************************/

#include <schnapps/core/map_handler.h>
#include <schnapps/core/schnapps.h>
#include <schnapps/core/view.h>

#include <cgogn/rendering/drawer.h>

namespace schnapps
{

MapHandlerGen::MapHandlerGen(const QString& name, SCHNApps* schnapps, MapBaseData* map) :
	name_(name),
	schnapps_(schnapps),
	map_(map),
	show_bb_(true),
	bb_diagonal_size_(.0f),
	bb_color_(Qt::green)
{
	connect(&frame_, SIGNAL(manipulated()), this, SLOT(frame_changed()));

	transformation_matrix_.setToIdentity();
}

MapHandlerGen::~MapHandlerGen()
{
	delete map_;
}

bool MapHandlerGen::is_selected_map() const
{
	return schnapps_->get_selected_map() == this;
}

/*********************************************************
 * MANAGE FRAME
 *********************************************************/

QMatrix4x4 MapHandlerGen::get_frame_matrix() const
{
	QMatrix4x4 m;
	GLdouble tmp[16];
	frame_.getMatrix(tmp);
	for (unsigned int i=0; i<4; ++i)
		for (unsigned int j=0; j<4; ++j)
			m(j,i) = tmp[i*4+j];
	return m;
}

void MapHandlerGen::frame_changed()
{
	emit(bb_changed());
}

/*********************************************************
 * MANAGE BOUNDING BOX
 *********************************************************/

void MapHandlerGen::set_show_bb(bool b)
{
	show_bb_ = b;
	foreach (View* view, views_)
		view->update();
}

void  MapHandlerGen::set_bb_color(const QString& color)
{
	bb_color_ = QColor(color);
	update_bb_drawer();
}

bool MapHandlerGen::get_transformed_bb(qoglviewer::Vec& bb_min, qoglviewer::Vec& bb_max)
{
	if (!bb_.is_initialized())
		return false;

	const VEC3& min = bb_.min();
	const VEC3& max = bb_.max();

	cgogn::geometry::BoundingBox<VEC3> bb;

	QVector4D v4 = transformation_matrix_ * QVector4D(min[0], min[1], min[2], 1.0f);
	qoglviewer::Vec v = qoglviewer::Vec(v4[0], v4[1], v4[2]);
	qoglviewer::Vec vt = frame_.inverseCoordinatesOf(v);
	bb.add_point(VEC3(vt[0], vt[1], vt[2]));

	v4 = transformation_matrix_ * QVector4D(max[0], min[1], min[2], 1.0f);
	v = qoglviewer::Vec(v4[0], v4[1], v4[2]);
	vt = frame_.inverseCoordinatesOf(v);
	bb.add_point(VEC3(vt[0], vt[1], vt[2]));

	v4 = transformation_matrix_ * QVector4D(min[0], max[1], min[2], 1.0f);
	v = qoglviewer::Vec(v4[0], v4[1], v4[2]);
	vt = frame_.inverseCoordinatesOf(v);
	bb.add_point(VEC3(vt[0], vt[1], vt[2]));

	v4 = transformation_matrix_ * QVector4D(min[0], min[1], max[2], 1.0f);
	v = qoglviewer::Vec(v4[0], v4[1], v4[2]);
	vt = frame_.inverseCoordinatesOf(v);
	bb.add_point(VEC3(vt[0], vt[1], vt[2]));

	v4 = transformation_matrix_ * QVector4D(max[0], max[1], min[2], 1.0f);
	v = qoglviewer::Vec(v4[0], v4[1], v4[2]);
	vt = frame_.inverseCoordinatesOf(v);
	bb.add_point(VEC3(vt[0], vt[1], vt[2]));

	v4 = transformation_matrix_ * QVector4D(max[0], min[1], max[2], 1.0f);
	v = qoglviewer::Vec(v4[0], v4[1], v4[2]);
	vt = frame_.inverseCoordinatesOf(v);
	bb.add_point(VEC3(vt[0], vt[1], vt[2]));

	v4 = transformation_matrix_ * QVector4D(min[0], max[1], max[2], 1.0f);
	v = qoglviewer::Vec(v4[0], v4[1], v4[2]);
	vt = frame_.inverseCoordinatesOf(v);
	bb.add_point(VEC3(vt[0], vt[1], vt[2]));

	v4 = transformation_matrix_ * QVector4D(max[0], max[1], max[2], 1.0f);
	v = qoglviewer::Vec(v4[0], v4[1], v4[2]);
	vt = frame_.inverseCoordinatesOf(v);
	bb.add_point(VEC3(vt[0], vt[1], vt[2]));

	bb_min = qoglviewer::Vec(bb.min()[0], bb.min()[1], bb.min()[2]);
	bb_max = qoglviewer::Vec(bb.max()[0], bb.max()[1], bb.max()[2]);

	return true;
}

void MapHandlerGen::draw_bb(View* view, const QMatrix4x4 &pm, const QMatrix4x4 &mm)
{
	bb_drawer_.call_list(pm, mm, view);
}

void MapHandlerGen::update_bb_drawer()
{
	if (bb_.is_initialized())
	{
		VEC3 bbmin = bb_.min();
		VEC3 bbmax = bb_.max();

		float shift = 0.005f * bb_diagonal_size_;
		bbmin -= VEC3(shift, shift, shift);
		bbmax += VEC3(shift, shift, shift);

		bb_drawer_.new_list();
		bb_drawer_.color3f(bb_color_.redF(), bb_color_.greenF(), bb_color_.blueF());
		bb_drawer_.line_width(2.0);
		bb_drawer_.begin(GL_LINE_LOOP);
			bb_drawer_.vertex3f(bbmin[0], bbmin[1], bbmin[2]);
			bb_drawer_.vertex3f(bbmin[0], bbmax[1], bbmin[2]);
			bb_drawer_.vertex3f(bbmax[0], bbmax[1], bbmin[2]);
			bb_drawer_.vertex3f(bbmax[0], bbmin[1], bbmin[2]);
		bb_drawer_.end();
		bb_drawer_.begin(GL_LINE_LOOP);
			bb_drawer_.vertex3f(bbmax[0], bbmax[1], bbmax[2]);
			bb_drawer_.vertex3f(bbmax[0], bbmin[1], bbmax[2]);
			bb_drawer_.vertex3f(bbmin[0], bbmin[1], bbmax[2]);
			bb_drawer_.vertex3f(bbmin[0], bbmax[1], bbmax[2]);
		bb_drawer_.end();
		bb_drawer_.begin(GL_LINES);
			bb_drawer_.vertex3f(bbmin[0], bbmin[1], bbmin[2]);
			bb_drawer_.vertex3f(bbmin[0], bbmin[1], bbmax[2]);
			bb_drawer_.vertex3f(bbmin[0], bbmax[1], bbmin[2]);
			bb_drawer_.vertex3f(bbmin[0], bbmax[1], bbmax[2]);
			bb_drawer_.vertex3f(bbmax[0], bbmax[1], bbmin[2]);
			bb_drawer_.vertex3f(bbmax[0], bbmax[1], bbmax[2]);
			bb_drawer_.vertex3f(bbmax[0], bbmin[1], bbmin[2]);
			bb_drawer_.vertex3f(bbmax[0], bbmin[1], bbmax[2]);
		bb_drawer_.end();
		bb_drawer_.end_list();
	}
}

/*********************************************************
 * MANAGE VBOs
 *********************************************************/

cgogn::rendering::VBO* MapHandlerGen::create_VBO(const QString& name)
{
	cgogn::rendering::VBO* vbo = get_VBO(name);
	if (!vbo)
	{
		vbo = new cgogn::rendering::VBO();
	}
}

cgogn::rendering::VBO* MapHandlerGen::get_VBO(const QString& name) const
{
	if (vbos_.contains(name))
		return vbos_[name];
	else
		return nullptr;
}

/*********************************************************
 * MANAGE LINKED VIEWS
 *********************************************************/

void MapHandlerGen::link_view(View* view)
{
	if (view && !views_.contains(view))
		views_.push_back(view);
}

void MapHandlerGen::unlink_view(View* view)
{
	views_.removeOne(view);
}

} // namespace schnapps
