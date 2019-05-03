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

#ifndef SCHNAPPS_PLUGIN_SURFACE_FEATURE_LINES_H_
#define SCHNAPPS_PLUGIN_SURFACE_FEATURE_LINES_H_

#include <schnapps/plugins/surface_feature_lines/dll.h>
#include <schnapps/plugins/surface_feature_lines/map_parameters.h>

#include <schnapps/core/types.h>
#include <schnapps/core/plugin_interaction.h>

#include <QAction>

namespace schnapps
{

class View;
namespace plugin_cmap_provider { class CMap2Handler; }

namespace plugin_surface_feature_lines
{


struct feature_point{
    cgogn::Dart d;
    double a;

    feature_point(){}
    feature_point(cgogn::Dart _d, double _a) : d(_d), a(_a){}

    bool is_valid(){
        if(!d.is_nil())
            return true;
        else
            return false;
    }
    static std::string cgogn_name_of_type() { return std::string("feature_line_struct"); }
    inline friend std::ostream& operator<<(std::ostream& out, const feature_point&)
    {
      return out;
    }
    inline friend std::istream& operator>> (std::istream& in, const feature_point&)
    {
      return in;
    }
};

using feature_line = std::vector<feature_point>;

class ComputeFeatureLines_Dialog;
using CMap2Handler = plugin_cmap_provider::CMap2Handler;

/**
* @brief Plugin for surface feature lines computation & display
*/
class SCHNAPPS_PLUGIN_SURFACE_FEATURE_LINES_API Plugin_SurfaceFeatureLines : public PluginInteraction
{
	Q_OBJECT
	Q_PLUGIN_METADATA(IID "SCHNApps.Plugin")
	Q_INTERFACES(schnapps::Plugin)

public:

	Plugin_SurfaceFeatureLines();
	inline ~Plugin_SurfaceFeatureLines() override {}
	static QString plugin_name();

	MapParameters& parameters(CMap2Handler* mh);

private:

	bool enable() override;
	void disable() override;

	inline void draw(View*, const QMatrix4x4&, const QMatrix4x4&) override {}
	void draw_object(View* view, Object* o, const QMatrix4x4& proj, const QMatrix4x4& mv) override;

	inline bool keyPress(View*, QKeyEvent*) override { return true; }
	inline bool keyRelease(View*, QKeyEvent*) override { return true; }
	inline bool mousePress(View*, QMouseEvent*) override { return true; }
	inline bool mouseRelease(View*, QMouseEvent*) override { return true; }
	inline bool mouseMove(View*, QMouseEvent*) override { return true; }
	inline bool wheelEvent(View*, QWheelEvent*) override { return true; }

	inline void resizeGL(View*, int, int) override {}

	void view_linked(View*) override;
	void view_unlinked(View*) override;

private slots:

	// slots called from SCHNApps signals
	void schnapps_closing();

	// slots called from View signals
	void object_linked(Object* o);
	void object_unlinked(Object* o);

	// slots called from action signals
	void open_compute_feature_lines_dialog();

private:

	void add_linked_map(View* view, CMap2Handler* mh);
	void remove_linked_map(View* view, CMap2Handler* mh);

private slots:

	// slots called from CMap2Handler signals
	void linked_map_attribute_changed(cgogn::Orbit orbit, const QString& name);

public:

	void compute_feature_lines(
		CMap2Handler* mh,
		const QString& position_attribute_name,
		const QString& Kmax_attribute_name,
		const QString& kmax_attribute_name,
		const QString& Kmin_attribute_name,
		const QString& kmin_attribute_name
    );

    std::vector<std::vector<feature_point>> feature_lines;

private:

	ComputeFeatureLines_Dialog* compute_feature_lines_dialog_;
	QAction* compute_feature_lines_action_;
	std::map<CMap2Handler*, MapParameters> parameter_set_;
};

} // namespace plugin_surface_feature_lines

} // namespace schnapps

#endif // SCHNAPPS_PLUGIN_SURFACE_FEATURE_LINES_H_
