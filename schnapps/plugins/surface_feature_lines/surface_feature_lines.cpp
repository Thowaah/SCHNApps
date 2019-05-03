
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

#include <schnapps/plugins/surface_selection/surface_selection.h>

namespace schnapps
{

namespace plugin_surface_feature_lines
{
using CMapCellsSetGen = plugin_cmap_provider::CMapCellsSetGen;

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

void add_lines(CMap2* map, CMap2::VertexAttribute<VEC3> position, CMap2::EdgeAttribute<feature_point>& edge_fp, std::vector<VEC3>& line_vec, CMap2::Vertex O,CMap2::Vertex v1, CMap2::Vertex v2, double scal1, double scal2){
    VEC3 line_start;
    VEC3 line_end;

    CMap2::Edge e1 = CMap2::Edge(O.dart);
    //CMap2::Edge e2 = CMap2::Edge(map->phi1(O.dart));
    CMap2::Edge e3 = CMap2::Edge(map->phi_1(O.dart));

    line_start = position[O]*(1-scal1) + (scal1)*position[v1];
    line_end = (1-scal2) * position[O] +  scal2*position[v2];

    edge_fp[e1] = feature_point(O.dart, scal1);
    edge_fp[e3] = feature_point(map->phi_1(O.dart), scal2);
    line_vec.push_back(line_start);
    line_vec.push_back(line_end);
}

//void singular_feature_line(CMap2* map, CMap2::VertexAttribute<VEC3>& K, CMap2::EdgeAttribute<feature_point>& edge_fp,CMap2::VertexAttribute<VEC3>& position,std::vector<VEC3>& lines){
//    map->foreach_cell([&] (CMap2::Face f){
//        int n_fl = 0;
//        std::vector<VEC3> fl_pts;
//        if(!isRegular(map,f,K)){
//           //affectation du point sur l'arete dans le cas singulier
//            // map->foreach_incident_edge(f,[&] (CMap2::Edge e){
//                //if(edge_fp[e].any()) n_fl++;
//                //fl_pts.push_back(edge_fp[e]);
//           // });
//            if(n_fl == 2){
//                lines.push_back(fl_pts[0]);
//                lines.push_back(fl_pts[1]);
//            } else if(n_fl == 3){
//                VEC3 barycentre;
//                map->foreach_incident_vertex(f, [&] (CMap2::Vertex v){
//                   //calcul barycentre
//                   barycentre = barycentre + position[v];
//                });
//                barycentre /= 3;

//                lines.push_back(fl_pts[0]);
//                lines.push_back(barycentre);
//                lines.push_back(fl_pts[1]);
//                lines.push_back(barycentre);
//                lines.push_back(fl_pts[2]);
//                lines.push_back(barycentre);
//            }
//        }
//    });
//}

void compute_feature_line(bool isKmax,CMap2* map, CMap2::VertexAttribute<VEC3>& K, CMap2::VertexAttribute<SCALAR>& e, CMap2::VertexAttribute<VEC3>& position, CMap2::VertexAttribute<SCALAR>& kmax, CMap2::VertexAttribute<SCALAR>& kmin, CMap2::FaceAttribute<SCALAR>& face_area, CMap2::VertexAttribute<SCALAR>& star_area, CMap2::FaceAttribute<VEC3>& grad_k, CMap2::EdgeAttribute<feature_point>& edge_fp,std::vector<VEC3>& lines){
    //vérifier si on a Kmin ou Kmax


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
//
            if( isKmax &&(grad_e.dot(sum_K)<0) && (std::abs(sum_kmax) > std::abs(sum_kmin))){
                if(e[v1] > 0 && e[v2] < 0 && e[v3] < 0){
                    add_lines(map,position, edge_fp,lines, v1, v2, v3, e[v1]/(e[v1] + -e[v2]) , e[v1]/(e[v1] + -e[v3]));
                } else if (e[v1] > 0 && e[v2] > 0 && e[v3] < 0){
                    add_lines(map, position,edge_fp,lines, v3, v1, v2, -e[v3]/(-e[v3] + e[v1]) , -e[v3]/(-e[v3] + e[v2]));
                } else if (e[v1] < 0 && e[v2] > 0 && e[v3] < 0){
                    add_lines(map, position,edge_fp,lines, v2, v1, v3, e[v2]/(e[v2] + (-e[v1])) , e[v2]/(e[v2] + (-e[v3])));
                } else if (e[v1] < 0 && e[v2] > 0 && e[v3] > 0){
                    add_lines(map, position,edge_fp,lines, v1, v2, v3, -e[v1]/(-e[v1] + e[v2]) , -e[v1]/(-e[v1] + e[v3]));
                } else if (e[v1] < 0 && e[v2] < 0 && e[v3] > 0){
                    add_lines(map, position,edge_fp,lines, v3, v1, v2, e[v3]/(e[v3] + -e[v1]) , e[v3]/(e[v3] + -e[v2]));
                } else  if(e[v1] > 0 && e[v2] < 0 && e[v3] > 0){
                    add_lines(map, position,edge_fp,lines, v2, v1, v3, -e[v2]/(-e[v2] + (e[v1])) , -e[v2]/(-e[v2] + (e[v3])));
                } else {
                    return;
                }
            } else if( !isKmax && (grad_e.dot(sum_K)>0) && (std::abs(sum_kmax) < std::abs(sum_kmin))){
                if(e[v1] > 0 && e[v2] < 0 && e[v3] < 0){
                    add_lines(map,position,edge_fp,lines, v1, v2, v3, e[v1]/(e[v1] + -e[v2]) , e[v1]/(e[v1] + -e[v3]));
                } else if (e[v1] > 0 && e[v2] > 0 && e[v3] < 0){
                    add_lines(map,position,edge_fp,lines, v3, v1, v2, -e[v3]/(-e[v3] + e[v1]) , -e[v3]/(-e[v3] + e[v2]));
                } else if (e[v1] < 0 && e[v2] > 0 && e[v3] < 0){
                    add_lines(map,position,edge_fp,lines, v2, v1, v3, e[v2]/(e[v2] + (-e[v1])) , e[v2]/(e[v2] + (-e[v3])));
                } else if (e[v1] < 0 && e[v2] > 0 && e[v3] > 0){
                    add_lines(map,position,edge_fp,lines, v1, v2, v3, -e[v1]/(-e[v1] + e[v2]) , -e[v1]/(-e[v1] + e[v3]));
                } else if (e[v1] < 0 && e[v2] < 0 && e[v3] > 0){
                    add_lines(map,position,edge_fp,lines, v3, v1, v2, e[v3]/(e[v3] + -e[v1]) , e[v3]/(e[v3] + -e[v2]));
                } else  if(e[v1] > 0 && e[v2] < 0 && e[v3] > 0){
                    add_lines(map,position,edge_fp,lines, v2, v1, v3, -e[v2]/(-e[v2] + (e[v1])) , -e[v2]/(-e[v2] + (e[v3])));
                } else {
                    return;
                }
            }
        }
    });
}

CMap2::Edge find_start(CMap2* map, CMap2::Edge _e, CMap2::EdgeAttribute<feature_point>& edge_fp){
    CMap2::CellMarker<CMap2::Edge::ORBIT> em(*map);
    em.mark(_e);
    CMap2::Edge first_valid_neighbour = _e;

    int n_valid = 0;

    do{

        //creer un marqueur
        n_valid = 0;



        map->foreach_incident_face(_e, [&] (CMap2::Face f){
            map->foreach_incident_edge(f, [&](CMap2::Edge e){
                if(edge_fp[e].is_valid() && !em.is_marked(e) && n_valid==0){
                     n_valid++;
                     first_valid_neighbour = e;
                     em.mark(e);
                }
            });


        });
    } while(n_valid !=0);

    return first_valid_neighbour;
}

void build_fl(CMap2* map, CMap2::Edge _e, feature_line& fl, CMap2::CellMarker<CMap2::Edge::ORBIT>& cm, CMap2::EdgeAttribute<feature_point>& edge_fp, CMapCellsSetGen* csg){
    CMap2::Edge next_e = find_start(map, _e, edge_fp);
    fl.push_back(edge_fp[next_e]);
    //map->foreach_incident_vertex(next_e, [&](CMap2::Vertex v){csg->select(v.dart);});
    cm.mark(next_e);

    std::vector<CMap2::Edge> marked_points;
    std::vector<CMap2::Edge> unmarked_points;
    //check if for every face there are no more unmarked edge
    int n_nope = 0;
    int max = 0;
    CMap2::Edge last_e;
    std::cout << "entering while loop in build_fl\n" << std::endl;
    do{
        last_e = next_e;
        max++;
        n_nope = 0;
        map->foreach_incident_face(next_e, [&] (CMap2::Face f){
            marked_points.clear();
            unmarked_points.clear();
            map->foreach_incident_edge(f, [&] (CMap2::Edge e){
                if(cm.is_marked(e) && edge_fp[e].is_valid()){
                    marked_points.push_back(e);
                } else if(!cm.is_marked(e) && edge_fp[e].is_valid()){
                    unmarked_points.push_back(e);
                }
                //on met juste les points contenant des fl, ceux représentant le côté restant n'y sont pas
            });
            if(marked_points.size()>1){
                n_nope++;
            } else if(marked_points.size()>0 && unmarked_points.size() > 0){
                for(CMap2::Edge e : unmarked_points){
                    fl.push_back(edge_fp[e]);
                    //map->foreach_incident_vertex(next_e, [&](CMap2::Vertex v){csg->select(v.dart);});
                    cm.mark(e);
                    next_e = e;
                }
            } else {
                n_nope++;
            }
        });
    }while(n_nope < 2 && max < 10000);
    if(max>=10000) std::cout << "!!!!!!!!!!!!!!!!while loop broken !!!!!!!!!!!!!!!!\n" << std::endl;
    std::cout << "end of build_fl" << std::endl;
}

void flood(CMap2* map, plugin_cmap_provider::CMapCellsSet<CMap2Handler, CMap2::Vertex>*  csg, CMapCellsSetGen* csg_outer, int spread){
    int n_current_vertex = 0;
    CMap2::CellMarker<CMap2::Vertex::ORBIT> cm(*map);
    std::deque<CMap2::Vertex> queue;
    csg->foreach_cell([&] (CMap2::Vertex v){
        queue.push_back(v);
        n_current_vertex++;
        cm.mark(v);
    });

    for (int i=0;i<spread;i++) {
        int n_vertex_next = 0;
        for(int i=0; i<n_current_vertex; i++){
            CMap2::Vertex& v_curr = queue.front();
            map->foreach_adjacent_vertex_through_edge(v_curr, [&] (CMap2::Vertex v){
                if(!cm.is_marked(v)){
                    csg_outer->select(v.dart);
                    cm.mark(v);
                    queue.push_back(v);
                    n_vertex_next++;
                }
            });
            queue.pop_front();
        }
        n_current_vertex = n_vertex_next;
    }
}

std::vector<feature_line> filter(const std::vector<feature_line> _feature_lines, unsigned long min_size){
    std::vector<feature_line> feature_lines(_feature_lines);
    for(unsigned long i=0; i<feature_lines.size(); i++){
        if(feature_lines[i].size()<min_size){
            feature_lines.erase(feature_lines.begin()+int(i));
        }
    }
    return feature_lines;
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
    CMapCellsSetGen* csg_fl = mh->add_cells_set(CMap2::Vertex::ORBIT,"feature_lines");
    csg_fl->set_mutually_exclusive(true);
    plugin_cmap_provider::CMapCellsSet<CMap2Handler, CMap2::Vertex> csg(*mh, "feature line");
    csg.set_mutually_exclusive(true);
    CMapCellsSetGen* csg_outer = mh->add_cells_set(CMap2::Vertex::ORBIT, "flooding");
    csg_outer->set_mutually_exclusive(true);
   // plugin_cmap_provider::CMapCellsSet<CMap2Handler, CMap2::Vertex> csgi_(*mh, "floodingtest");

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
    //feature point on edge
    //bug ici
    CMap2::EdgeAttribute<feature_point> edge_fp_min = map->add_attribute<feature_point, CMap2::Edge::ORBIT>("edge_fp_min");
    CMap2::EdgeAttribute<feature_point> edge_fp_max = map->add_attribute<feature_point, CMap2::Edge::ORBIT>("edge_fp_max");

    //pour afficher
    std::vector<VEC3> lines;

    std::cout << "Computing feature lines values\n" << std::endl;
    //appel à la fonction de calcul
    compute_feature_line(true,map,Kmax,emax,position,kmax,kmin,face_area,star_area,grad_kmax,edge_fp_max,lines);
   // feature_line(map,Kmax,emax,position,kmax,kmin,face_area,star_area,grad_kmax,edge_fl_max,lines);

    //calcul triangles singuliers
    //singular_feature_line(map,Kmax,edge_fp_max,position,lines);

    //construire les feature lines
    std::vector<feature_line> feature_lines_max;
    CMap2::CellMarker<CMap2::Edge::ORBIT> cm(*map);

    std::cout << "starting fl construction\n" << std::endl;
    //parcours de tous les sommets à la recherche d'un début de fl
    //et ajout des lignes detectées dans feature_lines_max
   // int stop = 0;
    map->foreach_cell([&] (CMap2::Face f){
        int n_fp = 0;

        //if(stop>0) return;

        CMap2::Edge em;
        CMap2::Edge e1(f.dart);
        CMap2::Edge e2(map->phi1(f.dart));
        CMap2::Edge e3(map->phi_1(f.dart));
        if(edge_fp_max[e1].is_valid() && !cm.is_marked(e1)){
            n_fp++;
            em=e1;
        }
        if(edge_fp_max[e2].is_valid() && !cm.is_marked(e2)){
            n_fp++;
            em=e2;
        }
        if(edge_fp_max[e3].is_valid() && !cm.is_marked(e3)){
            n_fp++;
            em=e3;
        }
//        map->foreach_incident_edge(f,[&] (CMap2::Edge e){
//            if(!edge_fp_max[e].d.is_nil()){
//                n_fp++;
//            }
//        });

        if(n_fp == 1){
            feature_line fl;

            //lancer parcours de la ligne ET marque les sommets
            std::cout << "building feature line\n" << std::endl;
            build_fl(map, em, fl, cm, edge_fp_max, &csg);

            //la mettre dans
            std::cout << "adding fl to vector\n" << std::endl;
            feature_lines_max.push_back(fl);
            //stop++;
        }
    });

    //enelever les fl plus petites que ...
    std::vector<feature_line> filtered_fls = filter(feature_lines_max, 50);

    for(auto& fl : filtered_fls){
        for(auto& fp : fl){
            map->foreach_incident_vertex(CMap2::Edge(fp.d), [&](CMap2::Vertex v){
                csg.select(v);
            });
        }
    }

    csg.foreach_cell([&] (CMap2::Vertex v){
       csg_fl->select(v.dart);
    });

    //inondation
    flood(map,&csg,csg_outer,5);




//    //affichage des lignes
    MapParameters& p = parameters(mh);
    p.update_lines(lines);

    for (View* view : mh->linked_views())
        view->update();
}

} // namespace plugin_surface_feature_lines

} // namespace schnapps
