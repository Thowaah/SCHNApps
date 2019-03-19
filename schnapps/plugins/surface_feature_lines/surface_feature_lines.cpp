
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

#include <cgogn/geometry/algos/area.h>
#include <cgogn/geometry/algos/normal.h>

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


//verifie la régularité d'un triangle
bool isRegular(CMap2* map, CMap2::Face f, CMap2::VertexAttribute<VEC3>& Ki)
{

        VEC3 v1 = Ki[CMap2::Vertex(f.dart)];
        VEC3 v2 = Ki[CMap2::Vertex(map->phi1(f.dart))];
        VEC3 v3 = Ki[CMap2::Vertex(map->phi_1(f.dart))];

        if (v1.dot(v2) > 0) {
            if(v1.dot(v3) > 0) {
                if(v2.dot(v3) > 0) {
                    return true;
                } else {
                    return false;
                }
            } else {
                if(v1.dot(-1*v3)>0){
                    if(v2.dot(-1*v3)>0){
                        //modifier Ki v3
                        Ki[map->phi_1(f.dart)] = -1*v3;
                        return true;
                    } else {
                        return false;
                    }
                } else {
                    return false;
                }
            }
        } else {
            if(v1.dot(-1*v2) > 0){
                if(v1.dot(v3) > 0) {
                    if((-1*v2).dot(v3) > 0) {
                        //modifier Ki v2
                        Ki[map->phi1(f.dart)] = -1*v2;
                        return true;
                    } else {
                        return false;
                    }
                } else {
                    if(v1.dot(-1*v3)>0){
                        if((-1*v2).dot(-1*v3)>0){
                            //modifier Ki v3 et v2
                            Ki[map->phi1(f.dart)] = -1*v2;
                            Ki[map->phi_1(f.dart)] = -1*v3;
                            return true;
                        } else {
                            return false;
                        }
                    } else {
                        return false;
                    }
                }
            } else {
                return false;
            }
        }
}

void add_lines(CMap2::VertexAttribute<VEC3> position, std::vector<VEC3>& line_vec, CMap2::Vertex O,CMap2::Vertex v1, CMap2::Vertex v2, double scal1, double scal2){
    VEC3 line_start;
    VEC3 line_end;

    line_start = position[O]*(1-scal1) + (scal1)*position[v1];
    line_end = (1-scal2) * position[O] +  scal2*position[v2];

    line_vec.push_back(line_start);
    line_vec.push_back(line_end);
}

void feature_line(CMap2* map, CMap2::VertexAttribute<VEC3>& K, CMap2::VertexAttribute<SCALAR>& e, CMap2::VertexAttribute<VEC3>& position, CMap2::VertexAttribute<SCALAR>& kmax, CMap2::VertexAttribute<SCALAR>& kmin, CMap2::FaceAttribute<SCALAR>& face_area, CMap2::VertexAttribute<SCALAR>& star_area, CMap2::FaceAttribute<VEC3>& grad_k, std::vector<VEC3>& lines){
    //calcul de la feature line sur chaque face sequentiellement
    map->foreach_cell([&] (CMap2::Face f){

        if(isRegular(map, f, K)){
            //pour les 3 vertex de la face
            map->foreach_incident_vertex(f, [&] (CMap2::Vertex v){

                //double emin = 0.0;
                double emax_val = 0.0;

                map->foreach_incident_face(v,[&] (CMap2::Face _f){
                    //emax += K[v].dot(grad_kmax[f]);
                    emax_val+= face_area[_f] * grad_k[_f].dot(K[v]);
                });

                e[v] = 1/star_area[v] * emax_val;
            });
            //fin calcul emax

            //vérification condition
            CMap2::Vertex v1 = CMap2::Vertex(f.dart);
            CMap2::Vertex v2 = CMap2::Vertex(map->phi1(f.dart));
            CMap2::Vertex v3 = CMap2::Vertex(map->phi_1(f.dart));

            VEC3 n = cgogn::geometry::normal(*map, f, position);

            VEC3 v12 = position[v2] - position[v1];
            VEC3 v23 = position[v3] - position[v2];
            VEC3 v31 = position[v1] - position[v3];

            VEC3 grad_e = n.cross(v12)/2*face_area[f] * e[v3] + n.cross(v23)/2*face_area[f] *  e[v1] + n.cross(v31)/2*face_area[f] * e[v2];

            VEC3 sum_K = VEC3(0.0,0.0,0.0);
            double sum_kmax = 0.0;
            double sum_kmin = 0.0;

            map->foreach_incident_vertex(f, [&] (CMap2::Vertex v){
                sum_K += K[v];
                sum_kmax += kmax[v];
                sum_kmin += kmin[v];
            });

            if((grad_e.dot(sum_K)<0) && (std::abs(sum_kmax) > std::abs(sum_kmin))){
                if(e[v1] > 0 && e[v2] < 0 && e[v3] < 0){
                    add_lines(position,lines, v1, v2, v3, e[v1]/(e[v1] + -e[v2]) , e[v1]/(e[v1] + -e[v3]));
                } else if (e[v1] > 0 && e[v2] > 0 && e[v3] < 0){
                    add_lines(position,lines, v3, v1, v2, -e[v3]/(-e[v3] + e[v1]) , -e[v3]/(-e[v3] + e[v2]));
                } else if (e[v1] < 0 && e[v2] > 0 && e[v3] < 0){
                    add_lines(position,lines, v2, v1, v3, e[v2]/(e[v2] + (-e[v1])) , e[v2]/(e[v2] + (-e[v3])));
                } else if (e[v1] < 0 && e[v2] > 0 && e[v3] > 0){
                    add_lines(position,lines, v1, v2, v3, -e[v1]/(-e[v1] + e[v2]) , -e[v1]/(-e[v1] + e[v3]));
                } else if (e[v1] < 0 && e[v2] < 0 && e[v3] > 0){
                    add_lines(position,lines, v3, v1, v2, e[v3]/(e[v3] + -e[v1]) , e[v3]/(e[v3] + -e[v2]));
                } else  if(e[v1] > 0 && e[v2] < 0 && e[v3] > 0){
                    add_lines(position,lines, v2, v1, v3, -e[v2]/(-e[v2] + (e[v1])) , -e[v2]/(-e[v2] + (e[v3])));
                } else {
                    return;
                }
            }
        }
    });
}


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


	// features pour les calculs
    CMap2::FaceAttribute<SCALAR> face_area = mh->add_attribute<SCALAR, CMap2::Face::ORBIT>("face_area");
    cgogn::geometry::compute_area<CMap2::Face>(*map, position, face_area);

	CMap2::VertexAttribute<SCALAR> star_area = mh->add_attribute<SCALAR, CMap2::Vertex::ORBIT>("star_area");
    cgogn::geometry::compute_area<CMap2::Vertex>(*map, position, star_area);

    //calcul du gradient sur chaque face
    CMap2::FaceAttribute<VEC3> grad_kmax = mh->add_attribute<VEC3, CMap2::Face::ORBIT>("grad_kmax");
    CMap2::FaceAttribute<VEC3> grad_kmin = mh->add_attribute<VEC3, CMap2::Face::ORBIT>("grad_kmin");


    map->foreach_cell([&] (CMap2::Face f)
    {
       CMap2::Vertex v1 = CMap2::Vertex(f.dart);
       CMap2::Vertex v2 = CMap2::Vertex(map->phi1(f.dart));
       CMap2::Vertex v3 = CMap2::Vertex(map->phi_1(f.dart));

       VEC3 n = cgogn::geometry::normal(*map, f, position);

       VEC3 v12 = position[v2] - position[v1];
       VEC3 v23 = position[v3] - position[v2];
       VEC3 v31 = position[v1] - position[v3];

       grad_kmax[f] = n.cross(v12)/2*face_area[f] * kmax[v3] + n.cross(v23)/2*face_area[f] *  kmax[v1] + n.cross(v31)/2*face_area[f] * kmax[v2];
       grad_kmin[f] = n.cross(v12)/2*face_area[f] * kmin[v3] + n.cross(v23)/2*face_area[f] *  kmin[v1] + n.cross(v31)/2*face_area[f] * kmin[v2];


    });

    //CMap2::FaceAttribute<VEC3> grad_kmin = mh->add_attribute<VEC3, CMap2::Face::ORBIT>("grad_kmin");

	//features à afficher
	CMap2::VertexAttribute<SCALAR> emin = mh->add_attribute<SCALAR, CMap2::Vertex::ORBIT>("emin");
	CMap2::VertexAttribute<SCALAR> emax = mh->add_attribute<SCALAR, CMap2::Vertex::ORBIT>("emax");

    //pour afficher
    std::vector<VEC3> lines;

    //appel à la fonction de calcul
    feature_line(map,Kmax,emax,position,kmax,kmin,face_area,star_area,grad_kmax,lines);


    //affichage des lignes
	MapParameters& p = parameters(mh);
	p.update_lines(lines);

	for (View* view : mh->linked_views())
		view->update();
}

} // namespace plugin_surface_feature_lines

} // namespace schnapps
