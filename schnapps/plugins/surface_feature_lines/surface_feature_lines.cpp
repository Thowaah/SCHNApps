
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
   // map->foreach_cell([&] (CMap2::Face f){
        //------------->récupération des sommets

        //VEC3 v1 = Kmax[CMap2::Vertex(f.dart)];
        //VEC3 v2 = Kmax[CMap2::Vertex(map->phi1(f.dart))];
        //VEC3 v3 = Kmax[CMap2::Vertex(map->phi_1(f.dart))];

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

        //isRegular[f] = true;
        //isRegular[f] = false;
    //});
}

void add_lines(CMap2::VertexAttribute<VEC3> position, std::vector<VEC3>& line_vec, CMap2::Vertex O,CMap2::Vertex v1, CMap2::Vertex v2, double scal1, double scal2){
    VEC3 line_start;
    VEC3 line_end;

    //line_start = position[O]*scal1 + (1-scal1)*position[v1];
    line_start = position[O]*(1-scal1) + (scal1)*position[v1];
    //line_end = scal2 * position[O] +  (1-scal2)*position[v2];
    line_end = (1-scal2) * position[O] +  scal2*position[v2];

    line_vec.push_back(line_start);
    line_vec.push_back(line_end);
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

    map->foreach_cell([&] (CMap2::Face f)
    {
       CMap2::Vertex v1 = CMap2::Vertex(f.dart);
       CMap2::Vertex v2 = CMap2::Vertex(map->phi1(f.dart));
       CMap2::Vertex v3 = CMap2::Vertex(map->phi_1(f.dart));

       VEC3 n = cgogn::geometry::normal(*map, f, position);

       VEC3 v12 = position[v2] - position[v1];
       VEC3 v23 = position[v3] - position[v2];
       VEC3 v31 = position[v1] - position[v3];

       //grad_kmax[f] = n.cross(v12)/2*face_area[f] * kmax[v1] + n.cross(v23)/2*face_area[f] *  kmax[v2];
       grad_kmax[f] = n.cross(v12)/2*face_area[f] * kmax[v3] + n.cross(v23)/2*face_area[f] *  kmax[v1] + n.cross(v31)/2*face_area[f] * kmax[v2];

    });

    CMap2::FaceAttribute<VEC3> grad_kmin = mh->add_attribute<VEC3, CMap2::Face::ORBIT>("grad_kmin");
//idem

	//features à afficher
	CMap2::VertexAttribute<SCALAR> emin = mh->add_attribute<SCALAR, CMap2::Vertex::ORBIT>("emin");
	CMap2::VertexAttribute<SCALAR> emax = mh->add_attribute<SCALAR, CMap2::Vertex::ORBIT>("emax");

    CMap2::FaceAttribute<bool> regularity = mh->add_attribute<bool, CMap2::Face::ORBIT>("regularity");

    //pour afficher
    std::vector<VEC3> lines;

    //calcul de la feature line sur chaque face sequentiellement
    map->foreach_cell([&] (CMap2::Face f){

        if(isRegular(map, f, Kmax)){
            regularity[f] = true;
            //pour les 3 vertex de la face
            map->foreach_incident_vertex(f, [&] (CMap2::Vertex v){

                //double emin = 0.0;
                double emax_val = 0.0;

                map->foreach_incident_face(v,[&] (CMap2::Face _f){
                    //emax += Kmax[v].dot(grad_kmax[f]);
                    emax_val+= face_area[_f] * grad_kmax[_f].dot(Kmax[v]);
                });

                emax[v] = 1/star_area[v] * emax_val;
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

            VEC3 grad_emax = n.cross(v12)/2*face_area[f] * emax[v3] + n.cross(v23)/2*face_area[f] *  emax[v1] + n.cross(v31)/2*face_area[f] * emax[v2];

            VEC3 sum_Kmax = VEC3(0.0,0.0,0.0);
            double sum_kmax = 0.0;
            double sum_kmin = 0.0;

            map->foreach_incident_vertex(f, [&] (CMap2::Vertex v){
                sum_Kmax += Kmax[v];
                sum_kmax += kmax[v];
                sum_kmin += kmin[v];
            });

//void add_lines(std::vector<VEC3>& line_vec, VEC3 O,VEC3 v1, VEC3 v2, double scal1, double scal2)
            if((grad_emax.dot(sum_Kmax)<0) && (std::abs(sum_kmax) > std::abs(sum_kmin))){
                if(emax[v1] > 0 && emax[v2] < 0 && emax[v3] < 0){
                    //add_lines(position,lines, v1, v2, v3, -emax[v1]/(-emax[v1] + emax[v2]) , 1-(-emax[v1]/(-emax[v1] + emax[v3])));
                    add_lines(position,lines, v1, v2, v3, emax[v1]/(emax[v1] + -emax[v2]) , emax[v1]/(emax[v1] + -emax[v3]));
                } else if (emax[v1] > 0 && emax[v2] > 0 && emax[v3] < 0){
                    add_lines(position,lines, v3, v1, v2, -emax[v3]/(-emax[v3] + emax[v1]) , -emax[v3]/(-emax[v3] + emax[v2]));
                } else if (emax[v1] < 0 && emax[v2] > 0 && emax[v3] < 0){
                    add_lines(position,lines, v2, v1, v3, emax[v2]/(emax[v2] + (-emax[v1])) , emax[v2]/(emax[v2] + (-emax[v3])));
                } else if (emax[v1] < 0 && emax[v2] > 0 && emax[v3] > 0){
                    add_lines(position,lines, v1, v2, v3, -emax[v1]/(-emax[v1] + emax[v2]) , -emax[v1]/(-emax[v1] + emax[v3]));
                } else if (emax[v1] < 0 && emax[v2] < 0 && emax[v3] > 0){
                    add_lines(position,lines, v3, v1, v2, emax[v3]/(emax[v3] + -emax[v1]) , emax[v3]/(emax[v3] + -emax[v2]));
                } else  if(emax[v1] > 0 && emax[v2] < 0 && emax[v3] > 0){
                    add_lines(position,lines, v2, v1, v3, -emax[v2]/(-emax[v2] + (emax[v1])) , -emax[v2]/(-emax[v2] + (emax[v3])));
                } else {
                    return;
                }
            }
        } else {
            regularity[f] = false;
             }
    });




//    //calcul e_i
//    map->foreach_cell([&] (CMap2::Vertex v){
//        //double emin = 0.0;
//        double emax_val = 0.0;

//        map->foreach_incident_face(v,[&] (CMap2::Face f){
//            //emax += Kmax[v].dot(grad_kmax[f]);
//            emax_val+= face_area[f] * grad_kmax[f].dot(Kmax[v]);
//        });

//        emax[v] = 1/star_area[v] * emax_val;
//        //
//        //-------------->récupération de l'aire dans les attributs
//        //double area;// = ...;

//        //----------->pour chaque triangle qui touche v
//            //-------------->récupération de l'aire dans les attributs
//            //-------------->récupération du gradient dans les attributs
//            //-------------->récupération du vecteur Ki dans les attributs



//            //e_i += aire_triangle * produiduit scalaire du gradient_i de T et vecteur K_i

//        //pour tout
//       // emin = 1/star_area[v] * emin;
//    });

    //parcours
	// map->foreach_cell([&] (CMap2::Face f){
	// 	VEC3 c(0,0,0);
	// 	//bla[f] = VEC3();
	// 	//f.dart[f] = VEC3();
	// 	map->foreach_incident_vertex(f, [&] (CMap2::Vertex v){
    // 		c+=position[v];gedit
	// 	});
	// 	//bla[f] = c;
	// });
	// //inverse
	// map->foreach_cell([&] (CMap2::Vertex v){
	// 	VEC3 c(0,0,0);
	// 	//bla[f] = VEC3();
	// 	//f.dart[f] = VEC3();
	// 	map->foreach_incident_vertex(v, [&] (CMap2::Face f){
	// 		c+=bla[f];
	// 	});
	// 	//bla[v] = c;
	// });


	//verification de la régularité de la face


	//scaling kmin/kmax for faces
    //map->foreach_cell([&] (CMap2::Vertex v){
        //calcul de l'aire du 1-voisinage ? cellule de voronoi ?
        //double area = 0.0;

            //--------->parcours des voisins ?

		//kmin et kmax piecewise
        //double kmin_pw = 3/area * kmin; //conversion ?
        //double kmax_pw = 3/area * kmax;
		//-------------->affecter ces valeurs à des attributs de sommets
    //});

	//calcul du gradient de k_i
    //map->foreach_cell([&] (CMap2::Face f){
        //récupération des kmin_pw/kmax_pw des sommets de la face

      //  VEC3 v1;
      //  VEC3 v2;
      //  VEC3 v3;

        //calcul aire
      //  VEC3 v12 = v2 - v1;
      //  VEC3 v23 = v3 - v2;
      //  VEC3 v31 = v1 - v3;

        //double d12 = std::sqrt(pow(v12[0],2) + pow(v12[1],2) + pow(v12[2],2));
        //double d23 = std::sqrt(pow(v23[0],2) + pow(v23[1],2) + pow(v23[2],2));
        //double d31 = std::sqrt(pow(v31[0],2) + pow(v31[1],2) + pow(v31[2],2));

        //double s = 1/2 * (d12+d23+d31);
        //double area = std::sqrt(s*(s-d12)*(s-d23)*(s-d31));

        //stocker l'aire ?

		//calcul gradient
		//---------> récupérer normale
        //VEC3 normal;
//        normal * (v1 - v3) / 2*area + normal * (v3 - v2) / 2*area + normal * (v2 - v1) / 2*area;

		//affecter la valeur du gradient à la face
//	});



	// compute the feature lines
    //std::vector<VEC3> lines;
    //lines.push_back(mh->bb().max());
    //lines.push_back(mh->bb().min());

	// update render data in corresponding MapParameters
	MapParameters& p = parameters(mh);
	p.update_lines(lines);

	for (View* view : mh->linked_views())
		view->update();
}

} // namespace plugin_surface_feature_lines

} // namespace schnapps
