#include <boost/python.hpp>
#include "raven.hpp"
#include "numpy_utils.hpp"

using namespace RavenBSP;
using namespace trajopt;
using namespace Eigen;
using namespace OpenRAVE;
using std::vector;

namespace py = boost::python;

bool gInteractive = false;
py::object openravepy;



EnvironmentBasePtr GetCppEnv(py::object py_env) {
	py::object openravepy = py::import("openravepy");
	int id = py::extract<int>(openravepy.attr("RaveGetEnvironmentId")(py_env));
	EnvironmentBasePtr cpp_env = RaveGetEnvironment(id);
	return cpp_env;
}
KinBodyPtr GetCppKinBody(py::object py_kb, EnvironmentBasePtr env) {
	int id = py::extract<int>(py_kb.attr("GetEnvironmentId")());
	return env->GetBodyFromEnvironmentId(id);
}
KinBody::LinkPtr GetCppLink(py::object py_link, EnvironmentBasePtr env) {
	KinBodyPtr parent = GetCppKinBody(py_link.attr("GetParent")(), env);
	int idx = py::extract<int>(py_link.attr("GetIndex")());
	return parent->GetLinks()[idx];
}

class PyGraphHandle {
	vector<GraphHandlePtr> m_handles;
public:
	PyGraphHandle(const vector<GraphHandlePtr>& handles) : m_handles(handles) {}
	PyGraphHandle(GraphHandlePtr handle) : m_handles(1, handle) {}
	void SetTransparency1(float alpha) {
		BOOST_FOREACH(GraphHandlePtr& handle, m_handles) {
			SetTransparency(handle, alpha);
		}
	}
};

class PyOSGViewer {
public:
	PyOSGViewer(OSGViewerPtr viewer) : m_viewer(viewer) {}
	int Step() {
		m_viewer->UpdateSceneData();
		m_viewer->Draw();
		return 0;
	}
	void UpdateSceneData() {
		m_viewer->UpdateSceneData();
	}
	PyGraphHandle PlotKinBody(py::object py_kb) {
		return PyGraphHandle(m_viewer->PlotKinBody(GetCppKinBody(py_kb, m_viewer->GetEnv())));
	}
	PyGraphHandle PlotLink(py::object py_link) {
		return PyGraphHandle(m_viewer->PlotLink(GetCppLink(py_link, m_viewer->GetEnv())));
	}
	void SetTransparency(py::object py_kb, float alpha) {
		m_viewer->SetTransparency(GetCppKinBody(py_kb, m_viewer->GetEnv()), alpha);
	}
	void SetAllTransparency(float a) {
		m_viewer->SetAllTransparency(a);
	}
	void Idle() {
		assert(!!m_viewer);
		m_viewer->Idle();
	}
	PyGraphHandle DrawText(std::string text, float x, float y, float fontsize, py::object pycolor) {
		OpenRAVE::Vector color = OpenRAVE::Vector(py::extract<float>(pycolor[0]), py::extract<float>(pycolor[1]), py::extract<float>(pycolor[2]), py::extract<float>(pycolor[3]));
		return PyGraphHandle(m_viewer->drawtext(text, x, y, fontsize, color));
	}

private:
	OSGViewerPtr m_viewer;
	PyOSGViewer() {}
};
PyOSGViewer PyGetViewer(py::object py_env) {
	EnvironmentBasePtr env = GetCppEnv(py_env);
	OSGViewerPtr viewer = OSGViewer::GetOrCreate(env);
	ALWAYS_ASSERT(!!viewer);
	return PyOSGViewer(viewer);
}

class PyRavenBSPWrapper {
public:
	RavenBSPWrapperPtr wrapper;
	PyRavenBSPWrapper(py::object env) : wrapper(new RavenBSPWrapper()) {
		wrapper->env = GetCppEnv(env);
	}

	void set_start(py::object vec);
	void set_start_sigma(py::object mat);
	void set_goal_trans(py::object tf);
	void set_T(int T);

	void set_sigma_pts_scale(double scale);
	void set_insertion_factor(double factor);

	void set_controls(py::list list_of_vec);
	py::list get_controls();

	void set_manip_name(const string& manip_name);
	void set_link_name(const string& link_name);

	void createViewer(bool sim_plotting=false, bool stage_plotting=false) {
		OSGViewerPtr viewer = OSGViewer::GetOrCreate(wrapper->env);
		ALWAYS_ASSERT(!!viewer);
		wrapper->setViewer(viewer, sim_plotting, stage_plotting);
	}
	PyOSGViewer getViewer() {
		return PyOSGViewer(wrapper->viewer);
	}

	void initialize();
	bool finished();
	void solve();
	void simulate_execution();

	void run();
};

void PyRavenBSPWrapper::set_start(py::object py_vec) {
	wrapper->start = toEigen<RavenBSP::StateT>(py_vec);
}
void PyRavenBSPWrapper::set_start_sigma(py::object py_mat) {
	wrapper->start_sigma = toEigen<RavenBSP::VarianceT>(py_mat);
}
void PyRavenBSPWrapper::set_controls(py::list list_of_vec) {
	wrapper->controls.clear();
	int n = py::len(list_of_vec);
	for (int vec_ind=0; vec_ind < n; vec_ind++) {
		py::object py_vec = py::extract<py::object>(list_of_vec[vec_ind]);
		wrapper->controls.push_back(toEigen<RavenBSP::ControlT>(py_vec));
	}
}

py::list PyRavenBSPWrapper::get_controls() {
	py::list l;
	for (int i=0;i<wrapper->controls.size();i++) {
		l.append(toNdarray1(wrapper->controls[i].data(),wrapper->controls[i].size()));
	}
	return l;
}

void PyRavenBSPWrapper::set_goal_trans(py::object tf) {
	wrapper->goal_trans = toEigen<Matrix4d>(tf);
}
void PyRavenBSPWrapper::set_T(int T) {
	wrapper->T = T;
}
void PyRavenBSPWrapper::set_sigma_pts_scale(double scale) {
	wrapper->sigma_pts_scale = scale;
}
void PyRavenBSPWrapper::set_insertion_factor(double factor) {
	wrapper->insertion_factor = factor;
}

void PyRavenBSPWrapper::set_manip_name(const string& manip_name) {
	wrapper->manip_name = manip_name;
}
void PyRavenBSPWrapper::set_link_name(const string& link_name) {
	wrapper->link_name = link_name;
}

void PyRavenBSPWrapper::initialize() {
	wrapper->initialize();
}
bool PyRavenBSPWrapper::finished() {
	return wrapper->finished();
}
void PyRavenBSPWrapper::solve() {
	wrapper->solve();
}
void PyRavenBSPWrapper::simulate_execution() {
	wrapper->simulate_execution();
}

void PyRavenBSPWrapper::run() {
	wrapper->run();
}

BOOST_PYTHON_MODULE(cravenbsppy) {
	np_mod = py::import("numpy");
	py::object openravepy = py::import("openravepy");
	string pyversion = py::extract<string>(openravepy.attr("__version__"));
	if (OPENRAVE_VERSION_STRING != pyversion) {
		PRINT_AND_THROW("the openrave on your pythonpath is different from the openrave version that trajopt links to!");
	}

	py::class_< PyRavenBSPWrapper >("RavenBSPWrapper", py::init<py::object>())
			.def("set_start", &PyRavenBSPWrapper::set_start)
			.def("set_start_sigma", &PyRavenBSPWrapper::set_start_sigma)
			.def("set_goal_trans", &PyRavenBSPWrapper::set_goal_trans)
			.def("set_T", &PyRavenBSPWrapper::set_T)
			.def("set_controls", &PyRavenBSPWrapper::set_controls)
			.def("get_controls", &PyRavenBSPWrapper::get_controls)
			.def("set_manip_name", &PyRavenBSPWrapper::set_manip_name)
			.def("set_link_name", &PyRavenBSPWrapper::set_link_name)
			.def("initialize", &PyRavenBSPWrapper::initialize)
			.def("finished", &PyRavenBSPWrapper::finished)
			.def("solve", &PyRavenBSPWrapper::solve)
			.def("simulate_execution", &PyRavenBSPWrapper::simulate_execution)
			.def("run", &PyRavenBSPWrapper::run)
			;


	py::class_< PyGraphHandle >("GraphHandle", py::no_init)
			   .def("SetTransparency", &PyGraphHandle::SetTransparency1)
			   ;

	py::class_< PyOSGViewer >("OSGViewer", py::no_init)
			   .def("UpdateSceneData", &PyOSGViewer::UpdateSceneData)
			   .def("Step", &PyOSGViewer::Step)
			   .def("PlotKinBody", &PyOSGViewer::PlotKinBody)
			   .def("PlotLink", &PyOSGViewer::PlotLink)
			   .def("SetTransparency", &PyOSGViewer::SetTransparency)
			   .def("SetAllTransparency", &PyOSGViewer::SetAllTransparency)
			   .def("Idle", &PyOSGViewer::Idle)
			   .def("DrawText", &PyOSGViewer::DrawText)
			   ;
	py::def("GetViewer", &PyGetViewer, "Get OSG viewer for environment or create a new one");
}
