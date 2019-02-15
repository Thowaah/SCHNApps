
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

#include <schnapps/plugins/surface_feature_lines/surface_feature_lines.h>
#include <schnapps/plugins/surface_feature_lines/dialog_compute_feature_lines.h>

#include <schnapps/core/schnapps.h>

namespace schnapps
{

namespace plugin_surface_feature_lines
{

Plugin_SurfaceFeatureLines::Plugin_SurfaceFeatureLines()
{
	this->name_ = SCHNAPPS_PLUGIN_NAME;
}

QString Plugin_SurfaceFeatureLines::plugin_name()
{
	return SCHNAPPS_PLUGIN_NAME;
}

MapParameters& Plugin_SurfaceFeatureLines::parameters(CMap2Handler* mh)
{
	cgogn_message_assert(mh, "Try to access parameters for null map handler");

	if (parameter_set_.count(mh) == 0)
	{
		MapParameters& p = parameter_set_[mh];
//		p.initialize_gl();
		return p;
	}
	else
		return parameter_set_[mh];
}

bool Plugin_SurfaceFeatureLines::enable()
{
	compute_feature_lines_dialog_ = new ComputeFeatureLines_Dialog(schnapps_, this);
	compute_feature_lines_action_ = schnapps_->add_menu_action("Surface;Compute Feature Lines", "compute feature lines");
	connect(compute_feature_lines_action_, SIGNAL(triggered()), this, SLOT(open_compute_feature_lines_dialog()));

	connect(schnapps_, SIGNAL(schnapps_closing()), this, SLOT(schnapps_closing()));

	//charger un pointeur sur le plugin d'import
	//voir dans dialog feature lines (schnapps->enable_plugin)
	return true;
}

void Plugin_SurfaceFeatureLines::disable()
{
	disconnect(schnapps_, SIGNAL(schnapps_closing()), this, SLOT(schnapps_closing()));

	disconnect(compute_feature_lines_action_, SIGNAL(triggered()), this, SLOT(open_compute_feature_lines_dialog()));
	schnapps_->remove_menu_action(compute_feature_lines_action_);
	delete compute_feature_lines_dialog_;
}

void Plugin_SurfaceFeatureLines::draw_object(View* view, Object* o, const QMatrix4x4& proj, const QMatrix4x4& mv)
{
	CMap2Handler* mh = qobject_cast<CMap2Handler*>(o);
	if (mh)
	{
		view->makeCurrent();
		QOpenGLFunctions* ogl = QOpenGLContext::currentContext()->functions();

		const MapParameters& p = parameters(mh);

		if (p.nb_feature_lines() > 0)
		{
			p.shader_bold_line_param_feature_lines_->bind(proj, mv);
			ogl->glDrawArrays(GL_LINES, 0, p.nb_feature_lines() * 2);
			p.shader_bold_line_param_feature_lines_->release();
		}
	}
}

void Plugin_SurfaceFeatureLines::view_linked(View* view)
{
	connect(view, SIGNAL(object_linked(Object*)), this, SLOT(object_linked(Object*)));
	connect(view, SIGNAL(object_unlinked(Object*)), this, SLOT(object_unlinked(Object*)));

	for (Object* o : view->linked_objects())
	{
		CMap2Handler* mh = qobject_cast<CMap2Handler*>(o);
		if (mh)
			add_linked_map(view, mh);
	}
}

void Plugin_SurfaceFeatureLines::view_unlinked(View* view)
{
	disconnect(view, SIGNAL(object_linked(Object*)), this, SLOT(object_linked(Object*)));
	disconnect(view, SIGNAL(object_unlinked(Object*)), this, SLOT(object_unlinked(Object*)));

	for (Object* o : view->linked_objects())
	{
		CMap2Handler* mh = qobject_cast<CMap2Handler*>(o);
		if (mh)
			remove_linked_map(view, mh);
	}
}

void Plugin_SurfaceFeatureLines::object_linked(Object* o)
{
	View* view = static_cast<View*>(sender());
	CMap2Handler* mh = qobject_cast<CMap2Handler*>(o);
	if (mh)
		add_linked_map(view, mh);
}

void Plugin_SurfaceFeatureLines::add_linked_map(View* view, CMap2Handler* mh)
{
	connect(mh, SIGNAL(attribute_changed(cgogn::Orbit, const QString&)), this, SLOT(linked_map_attribute_changed(cgogn::Orbit, const QString&)), Qt::UniqueConnection);
}

void Plugin_SurfaceFeatureLines::object_unlinked(Object* o)
{
	View* view = static_cast<View*>(sender());
	CMap2Handler* mh = qobject_cast<CMap2Handler*>(o);
	if (mh)
		remove_linked_map(view, mh);
}

void Plugin_SurfaceFeatureLines::remove_linked_map(View* view, CMap2Handler* mh)
{
	disconnect(mh, SIGNAL(attribute_changed(cgogn::Orbit, const QString&)), this, SLOT(linked_map_attribute_changed(cgogn::Orbit, const QString&)));
}

void Plugin_SurfaceFeatureLines::linked_map_attribute_changed(cgogn::Orbit orbit, const QString& name)
{
	CMap2Handler* mh = qobject_cast<CMap2Handler*>(sender());

	if (orbit == CMap2::Vertex::ORBIT)
	{
		if (parameter_set_.count(mh) > 0ul)
		{
			MapParameters& p = parameter_set_[mh];
			// if the modified attribute is the attribute used as position
			// when computing feature lines for this surface, recompute feature lines
		}

		for (View* view : mh->linked_views())
			view->update();
	}
}

void Plugin_SurfaceFeatureLines::open_compute_feature_lines_dialog()
{
	compute_feature_lines_dialog_->show();
}

void Plugin_SurfaceFeatureLines::schnapps_closing()
{
	compute_feature_lines_dialog_->close();
}

/******************************************************************************/
/*                             PUBLIC INTERFACE                               */
/******************************************************************************/

void Plugin_SurfaceFeatureLines::compute_feature_lines(
	CMap2Handler* mh,
	const QString& position_attribute_name,
	const QString& Kmax_attribute_name,
	const QString& kmax_attribute_name,
	const QString& Kmin_attribute_name,
	const QString& kmin_attribute_name
)
{
	CMap2* map = mh->map();

	CMap2::VertexAttribute<VEC3> position = map->get_attribute<VEC3, CMap2::Vertex::ORBIT>(position_attribute_name.toStdString());
	if (!position.is_valid())
		return;

	CMap2::VertexAttribute<VEC3> Kmax = map->get_attribute<VEC3, CMap2::Vertex::ORBIT>(Kmax_attribute_name.toStdString());
	if (!Kmax.is_valid())
		return;

	CMap2::VertexAttribute<SCALAR> kmax = map->get_attribute<SCALAR, CMap2::Vertex::ORBIT>(kmax_attribute_name.toStdString());
	if (!kmax.is_valid())
		return;

	CMap2::VertexAttribute<VEC3> Kmin = map->get_attribute<VEC3, CMap2::Vertex::ORBIT>(Kmin_attribute_name.toStdString());
	if (!Kmin.is_valid())
		return;

	CMap2::VertexAttribute<SCALAR> kmin = map->get_attribute<SCALAR, CMap2::Vertex::ORBIT>(kmin_attribute_name.toStdString());
	if (!kmin.is_valid())
		return;

		//CMap2::FaceAttribute<VEC3> bla = ...
		//CMap2::VertexAttribute<VEC3> vbla = mh->add_attribute<VEC3, CMap2::Vertex::ORBIT>("bla");

	//parcours
	// map->foreach_cell([&] (CMap2::Face f){
	// 	VEC3 c(0,0,0);
	// 	//bla[f] = VEC3();
	// 	//f.dart[f] = VEC3();
	// 	map->foreach_incident_vertex(f, [&] (CMap2::Vertex v){
	// 		c+=position[v];
	// 	});
	// 	//bla[f] = c;
	// });
	// //inerse
	// map->foreach_cell([&] (CMap2::Vertex v){
	// 	VEC3 c(0,0,0);
	// 	//bla[f] = VEC3();
	// 	//f.dart[f] = VEC3();
	// 	map->foreach_incident_vertex(v, [&] (CMap2::Face f){
	// 		c+=bla[f];
	// 	});
	// 	//bla[v] = c;
	// });

	// compute the feature lines
	std::vector<VEC3> lines;
	lines.push_back(mh->bb().max());
	lines.push_back(mh->bb().min());

	// update render data in corresponding MapParameters
	MapParameters& p = parameters(mh);
	p.update_lines(lines);

	for (View* view : mh->linked_views())
		view->update();
}

} // namespace plugin_surface_feature_lines

} // namespace schnapps
