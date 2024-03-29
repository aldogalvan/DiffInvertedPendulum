
//------------------------------------------------------------------------------
#include "chai3d.h"
#include <GLFW/glfw3.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

//------------------------------------------------------------------------------
using namespace Eigen;
using namespace chai3d;
using namespace std;
//------------------------------------------------------------------------------

MatrixXd pApax()

{
    MatrixXd ret(4,4);
    ret << 0 , 0 , 0 , 1,
            0 , 0 , -1 , 0 ,
            0 , 1 , 0 , 0 ,
            -1 , 0 , 0 , 0;
    return ret;
     }


MatrixXd pApay()

{
    MatrixXd ret(4,4);
    ret << 0 , 0 , 1 , 0,
            0 , 0 , 0 , 1 ,
            -1 , 0 , 0 , 0 ,
            0 , -1 , 0 , 0;
    return ret;
}


MatrixXd pApaz()
{
    MatrixXd ret(4,4);
    ret << 0 , -1 , 0 , 0,
            1 , 0 , 0 , 0 ,
            0 , 0 , 0 , 1 ,
            0 , 0 , -1 , 0;
    return ret;
}

MatrixXd pApas()
{
    MatrixXd ret = MatrixXd::Identity(4,4);
    return ret;
}

Eigen::Quaterniond eulerToQuaternion(double roll, double pitch, double yaw) {
    // Eigen uses the following order of rotations by default: ZYX
    Eigen::Quaterniond q = Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ())
                           * Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY())
                           * Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX());
    return q;
}
VectorXd q2vec(Quaterniond q)
{
    return Vector4d(q.w(),q.x(),q.y(),q.z());
}
Quaterniond operator*(double scalar, const Quaterniond& q) {
    return Quaterniond(scalar * q.w(), scalar * q.x(), scalar * q.y(), scalar * q.z());
}

// Quaternion addition
Quaterniond operator+(const Quaterniond& q1, const Quaterniond& q2) {
    return Quaterniond(q1.w() + q2.w(), q1.x() + q2.x(), q1.y() + q2.y(), q1.z() + q2.z());
}



// Overload the << operator for Eigen::Quaterniond
std::ostream& operator<<(std::ostream& os, const Eigen::Quaterniond& q) {
    os << q.w() << " " << q.x() << " " << q.y() << " " << q.z() ;
    return os;
}


Matrix3d vec2skew(const Vector3d& v) {
    Matrix3d V;
    V << 0, -v.z(), v.y(),
            v.z(), 0,  -v.x(),
            -v.y(), v.x(), 0;
    return V;
}

MatrixXd Q(Quaterniond q)
{
    MatrixXd Q = MatrixXd::Zero(4,3);
    Q(0,0) = 0.5*q.w();
    Q(0,1) = 0.5*q.z();
    Q(0,2) = -0.5*q.y();
    Q(1,0) = -0.5*q.z();
    Q(1,1) = 0.5*q.w();
    Q(1,2) = 0.5*q.x();
    Q(2,0) = 0.5*q.y();
    Q(2,1) = -0.5*q.x();
    Q(2,2) = 0.5*q.w();
    Q(3,0) = -0.5*q.x();
    Q(3,1) = -0.5*q.y();
    Q(3,2) = -0.5*q.z();
    return Q;
}
MatrixXd pQpx()
{
    MatrixXd ret = MatrixXd::Zero(4,3);
    ret(1,2) = 1;
    ret(2,1) = -1;
    ret(3,0) = -1;
    return 0.5*ret;
}

MatrixXd pQpy()
{
    MatrixXd ret = MatrixXd::Zero(4,3);
    ret(0,2) = -1;
    ret(2,0) = 1;
    ret(3,1) = -1;
    return 0.5*ret;
}

MatrixXd pQpz()
{
    MatrixXd ret = MatrixXd::Zero(4,3);
    ret(0,1) = 1;
    ret(1,0) = -1;
    ret(3,2) = -1;
    return 0.5*ret;
}

MatrixXd pQpw()
{
    MatrixXd ret = MatrixXd::Zero(4,3);
    ret(0,0) = 1;
    ret(1,1) = 1;
    ret(2,2) = 1;
    return 0.5*ret;
}

MatrixXd pRpx(Quaterniond& q)
{
    MatrixXd ret(3,3);
    ret(0,0) = q.x(); ret(0,1) = q.y(); ret(0,2) = q.z();
    ret(1,0) = q.y(); ret(1,1) = -q.x(); ret(1,2) = -q.w();
    ret(2,0) = q.z(); ret(2,1) = q.w(); ret(2,2) = -q.x();
    return 2*ret;
}

MatrixXd pRpy(Quaterniond& q)
{
    MatrixXd ret(3,3);
    ret(0,0) = -q.y(); ret(0,1) = q.x(); ret(0,2) = q.w();
    ret(1,0) = q.x(); ret(1,1) = q.y(); ret(1,2) = q.z();
    ret(2,0) = -q.w(); ret(2,1) = q.z(); ret(2,2) = -q.y();
    return 2*ret;
}

MatrixXd pRpz(Quaterniond& q)
{
    MatrixXd ret(3,3);
    ret(0,0) = -q.z(); ret(0,1) = -q.w(); ret(0,2) = q.x();
    ret(1,0) = q.w(); ret(1,1) = -q.z(); ret(1,2) = q.y();
    ret(2,0) = q.x(); ret(2,1) = q.y(); ret(2,2) = q.z();
    return 2*ret;
}

MatrixXd pRpw(Quaterniond& q)
{
    MatrixXd ret(3,3);
    ret(0,0) = q.w(); ret(0,1) = -q.z(); ret(0,2) = q.y();
    ret(1,0) = q.z(); ret(1,1) = q.w(); ret(1,2) = -q.x();
    ret(2,0) = -q.y(); ret(2,1) = q.x(); ret(2,2) = q.w();
    return 2*ret;
}


MatrixXd q2mat(Quaterniond& q)
{
    MatrixXd C(4,4);
    C.block<3,3>(0,0) = vec2skew(q.vec());
    return C;
}

struct Pendulum
{
    Pendulum(cGenericTool* a_device, cShapeSphere* a_tool = NULL, cShapeSphere* a_marker = NULL)
    {
        device = a_device;
        tool = a_tool;
        marker = a_marker;

        // create a sphere (cursor) to represent the haptic device
        pendulum = new cMesh();
        cCreateCylinder(pendulum,l
                        ,r);
        pendulum->m_material->setRed();

        // align this pendulum with the center of mass
        cVector3d com(0,0,0);
        for (int i = 0 ; i < pendulum->getNumVertices(); i++)
        {
            com += pendulum->m_vertices->getLocalPos(i);
        }
        com /= pendulum->getNumVertices();
        pendulum->offsetVertices(-com);

        J = InertiaTensorCylinder(m,r,l);\
        J.setIdentity(); J *= 1e-3;
        Jinv = J.inverse();

        cVector3d pos = device->getDeviceGlobalPos();
        x = pos.eigen() - xc;
    }

    ~Pendulum(){}


    static Matrix3d InertiaTensorCylinder(double mass, double radius, double height) {
        Matrix3d inertiaTensor;
        double I_xx_yy = 1.0 / 12.0 * mass * (3.0 * radius * radius + height * height);
        double I_zz = 0.5 * mass * radius * radius;

        inertiaTensor << I_xx_yy, 0, 0,
                0, I_xx_yy, 0,
                0, 0, I_zz;

        return inertiaTensor;
    }

    void ImplicitEuler(double dt)
    {
        cout << " --------------------------------- " << endl;

        Matrix3d M_ = Matrix3d::Identity()*m;
        Vector3d v = P/m;
        MatrixXd R = q.matrix();
        Vector3d w = R*Jinv*R.transpose()*L;
        MatrixXd I = MatrixXd::Identity(3,3);
        Quaterniond qc = Quaterniond::Identity();

        xh = device->getDeviceGlobalPos().eigen();
        xh = Vector3d(0.005,0.005,0.005);
        vh = device->getDeviceGlobalLinVel().eigen();
        qh = Quaterniond(device->getDeviceGlobalRot().eigen());
        wh = device->getDeviceGlobalAngVel().eigen();

        cout << "xh = " << xh.transpose() << endl;
        cout << "qh = " << qh.w() << " " << qh.x() << " " << qh.y() << " " << qh.z() << endl;
        cout << "x = " << x.transpose() << endl;
        cout << "q = " << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << endl;
        cout << "v = " << v.transpose() << endl;
        cout << "w = " << w.transpose() << endl;

        Quaterniond q_inv = q.inverse();
        Quaterniond deltaq = (qh*qc.inverse()*q_inv).normalized();
        MatrixXd C = Vector4d(deltaq.x(),deltaq.y(),deltaq.z(),deltaq.w())*RowVector4d(q_inv.x(),q_inv.y(),q_inv.z(),q_inv.w());
        Vector3d u_c = 2 * acos(deltaq.w()) * deltaq.vec();

        cout << "u_c = " << u_c.transpose() << endl;

        // coupling force
        Fc = Kc*(xh - (x + q*xc)) + Bc*(vh - (v + w.cross(xc)));
        //Fc(0) = 0;
        Tc = ((Vector3d)(R*xc)).cross(Fc) + Kth*u_c + Bth*(wh - w);

        cout << "Fc = "<< Fc.transpose() << endl;
        cout << "Tc1 = " << ((Vector3d)(R*xc)).cross(Fc).transpose() << endl;
        cout << "Tc2 = " <<  Kth*u_c.transpose()  << endl;
        cout << "Tc = " << Tc.transpose() << endl;
        // sum total forces
        Vector3d F = Fc + m*Vector3d(0,0,-0.0981);
        Vector3d T = Tc;

        // coupling jacobian
        MatrixXd J_c = MatrixXd::Zero(13,13);

        MatrixXd pRpx_ = pRpx(q); //! CHECK!
        MatrixXd pRpy_ = pRpy(q); //! CHECK!
        MatrixXd pRpz_ = pRpz(q); //! CHECK!
        MatrixXd pRpw_ = pRpw(q); //! CHECK!

        // Derivative of angular velocity wrt quaternion
        VectorXd pwpqx = (pRpx_*Jinv*R.transpose() + R*Jinv*pRpx_.transpose())*L; //! CHECK
        VectorXd pwpqy = (pRpy_*Jinv*R.transpose() + R*Jinv*pRpy_.transpose())*L; //! CHECK
        VectorXd pwpqz = (pRpz_*Jinv*R.transpose() + R*Jinv*pRpz_.transpose())*L; //! CHECK
        VectorXd pwpqw = (pRpw_*Jinv*R.transpose() + R*Jinv*pRpw_.transpose())*L; //! CHECK

        // derivatives of force with respect to quaternion
        VectorXd pFcpqx= -Kc*pRpx_*xc + Bc * vec2skew(xc)*pwpqx; //! CHECK!
        VectorXd pFcpqy= -Kc*pRpy_*xc + Bc * vec2skew(xc)*pwpqy; //! CHECK!
        VectorXd pFcpqz= -Kc*pRpz_*xc + Bc * vec2skew(xc)*pwpqz; //! CHECK!
        VectorXd pFcpqw= -Kc*pRpw_*xc + Bc * vec2skew(xc)*pwpqw; //! CHECK!

        MatrixXd pupq;
        if (abs(1-deltaq.w()) > 1e-6) {
            pupq = 2 * acos(deltaq.w()) * C.block<3, 4>(0, 0) -
                   2 / (sqrt(1 - deltaq.w() * deltaq.w())) * deltaq.vec() * C.block<1, 4>(3, 0); //! CHECK!
        }
        else {
            pupq = MatrixXd::Zero(3,4); //! CHECK!
        }
        //        cout << "pupq = \n" << pupq << endl;
        // Elements of Jacobian
        MatrixXd pFcpx = -Kc*MatrixXd::Identity(3,3); //! CHECK!
        MatrixXd pTcpx = -Kc*vec2skew(R*xc); //! CHECK!
        MatrixXd pFcpq(3,4);
        pFcpq.col(0) = pFcpqx; //! CHECK!
        pFcpq.col(1) = pFcpqy; //! CHECK!
        pFcpq.col(2) = pFcpqz; //! CHECK!
        pFcpq.col(3) = pFcpqw; //! CHECK!
        MatrixXd pTcpq(3,4);
        pTcpq.col(0) = vec2skew(R*xc)*pFcpqx - vec2skew(Fc)*pRpx_*xc + Kth*pupq.col(0) - Bth*pwpqx; //! Check!
        pTcpq.col(1) = vec2skew(R*xc)*pFcpqy - vec2skew(Fc)*pRpy_*xc + Kth*pupq.col(1) - Bth*pwpqy; //! Check!
        pTcpq.col(2) = vec2skew(R*xc)*pFcpqz - vec2skew(Fc)*pRpz_*xc + Kth*pupq.col(2) - Bth*pwpqz; //! Check!
        pTcpq.col(3) = vec2skew(R*xc)*pFcpqw - vec2skew(Fc)*pRpw_*xc + Kth*pupq.col(3) - Bth*pwpqw; //! Check!
        MatrixXd pFcpP = -Bc*I/m; //! Check!
        MatrixXd pTcpP = -(Bc/m)*vec2skew(R*xc); //! Check!
        MatrixXd pFcpL = Bc*vec2skew(xc)*R*Jinv*R.transpose();
        MatrixXd pTcpL = (Bc*vec2skew(R*xc)*vec2skew(xc) - Bth*I)*R*Jinv*R.transpose();

        J_c.block<3,3>(7,0) = pFcpx;
        J_c.block<3,3>(10,0) = pTcpx;
        J_c.block<3,4>(7,3) = pFcpq;
        J_c.block<3,4>(10,3) = pTcpq;
        J_c.block<3,3>(7,7) = pFcpP;
        J_c.block<3,3>(10,7) = pTcpP;;
        J_c.block<3,3>(7,10) = pFcpL;
        J_c.block<3,3>(10,10) = pTcpL;

        MatrixXd Q_ = Q(q);
        MatrixXd pQpx_ = pQpx(); //! CHECK
        MatrixXd pQpy_ = pQpy(); //! CHECK
        MatrixXd pQpz_ = pQpz(); //! CHECK
        MatrixXd pQpw_ = pQpw(); //! CHECK

        MatrixXd pqdotpL = Q_*R*Jinv*R.transpose();
        MatrixXd pqdotpq(4,4);
        pqdotpq.col(0) = pQpx_*w + Q_*pwpqx;
        pqdotpq.col(1) = pQpy_*w + Q_*pwpqy;
        pqdotpq.col(2) = pQpz_*w + Q_*pwpqz;
        pqdotpq.col(3) = pQpw_*w + Q_*pwpqw;

        MatrixXd J_r = MatrixXd::Zero(13,13);
        J_r.block<3,3>(0,7) = I / m;
        J_r.block<4,4>(3,3) = pqdotpq;
        J_r.block<4,3>(3,10) = pqdotpL;

        // solve the system
        MatrixXd A = (MatrixXd::Identity(13,13) - dt*(J_c + J_r));
        VectorXd b = VectorXd::Zero(13);
        b.block<3,1>(0,0)  = dt * P / m;
        b.block<4,1>(3,0)  = dt * Q_ * w;
        b.block<3,1>(7,0)  = dt * F;
        b.block<3,1>(10,0) = dt * T;
        VectorXd y = A.fullPivLu().solve(b);
//        cout << "A = \n" << A << endl;
//        cout << "b = " << b.transpose() << endl;
        cout << "ds = " << y.transpose() << endl;

        // Implicit Euler Integration
        P += y.block<3,1>(7,0);
        L += y.block<3,1>(10,0);
        x += y.block<3,1>(0,0);
        //Quaterniond dq =  (dt/2)*Quaterniond(0,w(0),w(1),w(2))*Quaterniond(y(6),y(3),y(4),y(5));
        Quaterniond dq = Quaterniond(y(6),y(3),y(4),y(5));
        q = q + dq;
        cout << "dq = " << dq.w() << " " << dq.x() << " " << dq.y() << " " << dq.z() << endl;
        q.normalize();


        updateGraphics();
    }

    void ExplicitEuler(double dt)
    {
        cout << " --------------------------------- " << endl;

        Matrix3d M_ = Matrix3d::Identity()*m;
        Vector3d v = P/m;
        MatrixXd R = q.matrix();
        Vector3d w = R*Jinv*R.transpose()*L;
        MatrixXd I = MatrixXd::Identity(3,3);
        Quaterniond qc = Quaterniond::Identity();

        xh = device->getDeviceGlobalPos().eigen();
        vh = device->getDeviceGlobalLinVel().eigen();
        qh = Quaterniond(device->getDeviceGlobalRot().eigen());
        wh = device->getDeviceGlobalAngVel().eigen();

        cout << "xh = " << xh.transpose() << endl;
        cout << "x = " << x.transpose() << endl;
        cout << "q = " << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << endl;
        cout << "v = " << v.transpose() << endl;
        cout << "w = " << w.transpose() << endl;

        Quaterniond q_inv = q.inverse();
        Quaterniond deltaq = (qh*qc.inverse()*q_inv).normalized();
        Vector3d u_c = 2 * acos(deltaq.w()) * deltaq.vec();

        cout << "u_c = " << u_c.transpose() << endl;

        // coupling force
        Fc = Kc*(xh - (x + q*xc)) + Bc*(vh - (v + w.cross(xc)));
        Tc = ((Vector3d)(R*xc)).cross(Fc) + Kth*u_c + Bth*(wh - w);

        cout << "Fc  = "  << Fc.transpose() << endl;
        cout << "Tc1 = " << ((Vector3d)(R*xc)).cross(Fc).transpose() << endl;
        cout << "Tc2 = " <<  Kth*u_c.transpose()  << endl;
        cout << "Tc  = " << Tc.transpose() << endl;

        // Compute net force and torque
        Vector3d F = Fc + m*Vector3d(0,0,-0.0981);
        Vector3d T = Tc;

        // Update linear and angular momentum
        P += dt * F;
        L += dt * T;

        // Update position and orientation
        x += dt * P / m;
        w = R*Jinv*R.transpose()*L;
        Quaterniond qdot = dt * 0.5 * Quaterniond(0, w.x(), w.y(), w.z()) * q;
        cout << "qdot = " << qdot << endl;
        q.w() += qdot.w();
        q.x() += qdot.x();
        q.y() += qdot.y();
        q.z() += qdot.z();
        q.normalize();

        cVector3d temp = pendulum->getLocalPos();
        marker->setLocalPos(temp);

        updateGraphics();
    }

    void updateGraphics()
    {
        tool->setLocalPos(xh);
        pendulum->setLocalPos(x);
        pendulum->setLocalRot(q.matrix());
        cVector3d temp = pendulum->getLocalPos();
        marker->setLocalPos(temp);
    }

    // parameters for the pendulum
    double r = 0.01; // the radius (just for visualization)
    double l = 0.2; // the length
    double m = 0.1; // the mass
    MatrixXd J; // the moment of inertia
    MatrixXd Jinv; // the inverse moment of inertia

    // the tool
    cGenericTool* device;
    cShapeSphere* tool;
    Vector3d xc = Vector3d(0.,0.,-l/2); // the coupling position defined wrt pendulum
    Quaterniond qc = Quaterniond::Identity(); // the quaternion orientation defined wrt pendulum
    double Kc = 100; // coupling stiffness
    double Kth = 0; // angular coupling stiffness
    double Bc = 1; // linear damping
    double Bth = 0; // angular damping
    Vector3d Fc = Vector3d::Zero();
    Vector3d Tc = Vector3d::Zero();

    // the pendulum
    cMesh* pendulum;

    // pendulum states
    Quaterniond q = Quaterniond(AngleAxisd(45, Eigen::Vector3d::UnitX())); //Quaterniond::Identity(); // the orientation
    Vector3d L = Vector3d(0,0,0); // the angular velocity
    Vector3d x = Vector3d::Zero(); // the linear displacement
    Vector3d P = Vector3d::Zero(); // the linear velocity

    // haptic states
    Quaterniond qh = Quaterniond::Identity(); // the orientation
    Vector3d wh = Vector3d(0,0,0); // the angular velocity
    Vector3d xh = Vector3d::Zero(); // the linear displacement
    Vector3d vh = Vector3d::Zero(); // the linear velocity

    // a marker for debugging
    cShapeSphere* marker;
};

Pendulum* pendulum;

//------------------------------------------------------------------------------
// GENERAL SETTINGS
//------------------------------------------------------------------------------

cStereoMode stereoMode = C_STEREO_DISABLED;

// fullscreen mode
bool fullscreen = false;

// mirrored display
bool mirroredDisplay = false;

//------------------------------------------------------------------------------
// DECLARED VARIABLES
//------------------------------------------------------------------------------

// a world that contains all objects of the virtual environment
cWorld* world;

// a camera to render the world in the window display
cCamera* camera;

// a light source to illuminate the objects in the world
cDirectionalLight *light;

// a haptic device handler
cHapticDeviceHandler* handler;

// a pointer to the current haptic device
cGenericHapticDevicePtr hapticDevice;

// the cursor representing the haptic device
cGenericTool* device;
cShapeSphere* tool;
cShapeSphere* marker;

// a label to display the haptic device model
cLabel* labelHapticDeviceModel;

// a label to display the position [m] of the haptic device
cLabel* labelHapticDevicePosition;

// a font for rendering text
cFontPtr font;

// a label to display the rate [Hz] at which the simulation is running
cLabel* labelRates;

// a flag to indicate if the haptic simulation currently running
bool simulationRunning = false;

// a flag to indicate if the haptic simulation has terminated
bool simulationFinished = true;

// a frequency counter to measure the simulation graphic rate
cFrequencyCounter freqCounterGraphics;

// a frequency counter to measure the simulation haptic rate
cFrequencyCounter freqCounterHaptics;

// haptic thread
cThread* hapticsThread;

// a handle to window display context
GLFWwindow* window = NULL;

// current width of window
int width  = 0;

// current height of window
int height = 0;

// swap interval for the display context (vertical synchronization)
int swapInterval = 1;


//------------------------------------------------------------------------------
// DECLARED FUNCTIONS
//------------------------------------------------------------------------------

// callback when the window display is resized
void windowSizeCallback(GLFWwindow* a_window, int a_width, int a_height);

// callback when an error GLFW occurs
void errorCallback(int error, const char* a_description);

// callback when a key is pressed
void keyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods);

// this function renders the scene
void updateGraphics(void);

// this function contains the main haptics simulation loop
void updateHaptics(void);

// this function closes the application
void close(void);

// this function writes the data
void writeData(std::ofstream&, double);

//==============================================================================
/*
    DEMO:   01-mydevice.cpp

    This application illustrates how to program forces, torques and gripper
    forces to your haptic device.

    In this example the application opens an OpenGL window and displays a
    3D cursor for the device connected to your computer. If the user presses
    onto the user button (if available on your haptic device), the color of
    the cursor changes from blue to green.

    In the main haptics loop function  "updateHaptics()" , the position,
    orientation and user switch status are read at each haptic cycle.
    Force and torque vectors are computed and sent back to the haptic device.
*/
//==============================================================================

int main(int argc, char* argv[])
{
    //--------------------------------------------------------------------------
    // INITIALIZATION
    //--------------------------------------------------------------------------

    cout << endl;
    cout << "-----------------------------------" << endl;
    cout << "CHAI3D" << endl;
    cout << "Demo: 01-mydevice" << endl;
    cout << "Copyright 2003-2016" << endl;
    cout << "-----------------------------------" << endl << endl << endl;
    cout << "Keyboard Options:" << endl << endl;
    cout << "[1] - Enable/Disable potential field" << endl;
    cout << "[2] - Enable/Disable damping" << endl;
    cout << "[f] - Enable/Disable full screen mode" << endl;
    cout << "[m] - Enable/Disable vertical mirroring" << endl;
    cout << "[q] - Exit application" << endl;
    cout << endl << endl;


    //--------------------------------------------------------------------------
    // OPENGL - WINDOW DISPLAY
    //--------------------------------------------------------------------------

    // initialize GLFW library
    if (!glfwInit())
    {
        cout << "failed initialization" << endl;
        cSleepMs(1000);
        return 1;
    }

    // set error callback
    glfwSetErrorCallback(errorCallback);

    // compute desired size of window
    const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    int w = 0.8 * mode->height;
    int h = 0.5 * mode->height;
    int x = 0.5 * (mode->width - w);
    int y = 0.5 * (mode->height - h);

    // set OpenGL version
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

    // set active stereo mode
    if (stereoMode == C_STEREO_ACTIVE)
    {
        glfwWindowHint(GLFW_STEREO, GL_TRUE);
    }
    else
    {
        glfwWindowHint(GLFW_STEREO, GL_FALSE);
    }

    // create display context
    window = glfwCreateWindow(w, h, "CHAI3D", NULL, NULL);
    if (!window)
    {
        cout << "failed to create window" << endl;
        cSleepMs(1000);
        glfwTerminate();
        return 1;
    }

    // get width and height of window
    glfwGetWindowSize(window, &width, &height);

    // set position of window
    glfwSetWindowPos(window, x, y);

    // set key callback
    glfwSetKeyCallback(window, keyCallback);

    // set resize callback
    glfwSetWindowSizeCallback(window, windowSizeCallback);

    // set current display context
    glfwMakeContextCurrent(window);

    // sets the swap interval for the current display context
    glfwSwapInterval(swapInterval);

#ifdef GLEW_VERSION
    // initialize GLEW library
    if (glewInit() != GLEW_OK)
    {
        cout << "failed to initialize GLEW library" << endl;
        glfwTerminate();
        return 1;
    }
#endif


    //--------------------------------------------------------------------------
    // WORLD - CAMERA - LIGHTING
    //--------------------------------------------------------------------------

    // create a new world.
    world = new cWorld();

    // set the background color of the environment
    world->m_backgroundColor.setBlack();

    // create a camera and insert it into the virtual world
    camera = new cCamera(world);
    world->addChild(camera);

    // position and orient the camera
    camera->set( cVector3d (0.5, 0.0, 0.0),    // camera position (eye)
                 cVector3d (0.0, 0.0, 0.0),    // look at position (target)
                 cVector3d (0.0, 0.0, 1.0));   // direction of the (up) vector

    // set the near and far clipping planes of the camera
    camera->setClippingPlanes(0.01, 10.0);

    // set stereo mode
    camera->setStereoMode(stereoMode);

    // set stereo eye separation and focal length (applies only if stereo is enabled)
    camera->setStereoEyeSeparation(0.01);
    camera->setStereoFocalLength(0.5);

    // set vertical mirrored display mode
    camera->setMirrorVertical(mirroredDisplay);

    // create a directional light source
    light = new cDirectionalLight(world);

    // insert light source inside world
    world->addChild(light);

    // enable light source
    light->setEnabled(true);

    // define direction of light beam
    light->setDir(-1.0, 0.0, 0.0);

    //--------------------------------------------------------------------------
    // HAPTIC DEVICE
    //--------------------------------------------------------------------------

    // create a haptic device handler
    handler = new cHapticDeviceHandler();

    // get a handle to the first haptic device
    handler->getDevice(hapticDevice, 0);

    // calibrate device (if necessary)
    hapticDevice->calibrate();

    // retrieve information about the current haptic device
    cHapticDeviceInfo info = hapticDevice->getSpecifications();

    // if the device has a gripper, enable the gripper to simulate a user switch
    hapticDevice->setEnableGripperUserSwitch(true);

    // create a tool (cursor) and insert into the world
    device = new cGenericTool(world);
    world->addChild(device);
    tool = new cShapeSphere(0.015);
    world->addChild(tool);
    tool->m_material->setBlueLight();
    marker = new cShapeSphere(0.015);
    world->addChild(marker);
    marker->m_material->setYellow();

    // connect the haptic device to the virtual tool
    device->setHapticDevice(hapticDevice);

    // map the physical workspace of the haptic device to a larger virtual workspace.
    device->setWorkspaceRadius(1.0);

    // haptic forces are enabled only if small forces are first sent to the device;
    // this mode avoids the force spike that occurs when the application starts when
    // the tool is located inside an object for instance.
    device->setWaitForSmallForce(true);

    // start the haptic tool
    device->start();

    //--------------------------------------------------------------------------
    // WIDGETS
    //--------------------------------------------------------------------------

    // create a font
    font = NEW_CFONTCALIBRI20();

    // create a label to display the haptic device model
    labelHapticDeviceModel = new cLabel(font);
    camera->m_frontLayer->addChild(labelHapticDeviceModel);
    labelHapticDeviceModel->setText(info.m_modelName);

    // create a label to display the position of haptic device
    labelHapticDevicePosition = new cLabel(font);
    camera->m_frontLayer->addChild(labelHapticDevicePosition);

    // create a label to display the haptic and graphic rate of the simulation
    labelRates = new cLabel(font);
    camera->m_frontLayer->addChild(labelRates);


    //--------------------------------------------------------------------------
    // START SIMULATION
    //--------------------------------------------------------------------------

    pendulum = new Pendulum(device,tool,marker);
    world->addChild(pendulum->pendulum);
    // create a thread which starts the main haptics rendering loop
    hapticsThread = new cThread();
    hapticsThread->start(updateHaptics, CTHREAD_PRIORITY_HAPTICS);

    // setup callback when application exits
    atexit(close);


    //--------------------------------------------------------------------------
    // MAIN GRAPHIC LOOP
    //--------------------------------------------------------------------------

    // call window size callback at initialization
    windowSizeCallback(window, width, height);

    // main graphic loop
    while (!glfwWindowShouldClose(window))
    {
        // get width and height of window
        glfwGetWindowSize(window, &width, &height);

        // render graphics
        updateGraphics();

        // swap buffers
        glfwSwapBuffers(window);

        // process events
        glfwPollEvents();

        // signal frequency counter
        freqCounterGraphics.signal(1);
    }

    // close window
    glfwDestroyWindow(window);

    // terminate GLFW library
    glfwTerminate();

    // exit
    return 0;
}

//------------------------------------------------------------------------------

void windowSizeCallback(GLFWwindow* a_window, int a_width, int a_height)
{
    // update window size
    width  = a_width;
    height = a_height;

    // update position of label
    labelHapticDeviceModel->setLocalPos(20, height - 40, 0);

    // update position of label
    labelHapticDevicePosition->setLocalPos(20, height - 60, 0);
}

//------------------------------------------------------------------------------

void errorCallback(int a_error, const char* a_description)
{
    cout << "Error: " << a_description << endl;
}

//------------------------------------------------------------------------------

void keyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods)
{
    // filter calls that only include a key press
    if ((a_action != GLFW_PRESS) && (a_action != GLFW_REPEAT))
    {
        return;
    }

        // option - exit
    else if ((a_key == GLFW_KEY_ESCAPE) || (a_key == GLFW_KEY_Q))
    {
        glfwSetWindowShouldClose(a_window, GLFW_TRUE);
    }

        // option - toggle fullscreen
    else if (a_key == GLFW_KEY_F)
    {
        // toggle state variable
        fullscreen = !fullscreen;

        // get handle to monitor
        GLFWmonitor* monitor = glfwGetPrimaryMonitor();

        // get information about monitor
        const GLFWvidmode* mode = glfwGetVideoMode(monitor);

        // set fullscreen or window mode
        if (fullscreen)
        {
            glfwSetWindowMonitor(window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
            glfwSwapInterval(swapInterval);
        }
        else
        {
            int w = 0.8 * mode->height;
            int h = 0.5 * mode->height;
            int x = 0.5 * (mode->width - w);
            int y = 0.5 * (mode->height - h);
            glfwSetWindowMonitor(window, NULL, x, y, w, h, mode->refreshRate);
            glfwSwapInterval(swapInterval);
        }
    }

        // option - toggle vertical mirroring
    else if (a_key == GLFW_KEY_M)
    {
        mirroredDisplay = !mirroredDisplay;
        camera->setMirrorVertical(mirroredDisplay);
    }
}

//------------------------------------------------------------------------------

void close(void)
{
    // stop the simulation
    simulationRunning = false;

    // wait for graphics and haptics loops to terminate
    while (!simulationFinished) { cSleepMs(100); }

    // close haptic device
    hapticDevice->close();

    // delete resources
    delete hapticsThread;
    delete world;
    delete handler;
}

//------------------------------------------------------------------------------

void updateGraphics(void)
{
    /////////////////////////////////////////////////////////////////////
    // UPDATE WIDGETS
    /////////////////////////////////////////////////////////////////////


    // update haptic and graphic rate data
    labelRates->setText(cStr(freqCounterGraphics.getFrequency(), 0) + " Hz / " +
                        cStr(freqCounterHaptics.getFrequency(), 0) + " Hz");

    // update position of label
    labelRates->setLocalPos((int)(0.5 * (width - labelRates->getWidth())), 15);


    /////////////////////////////////////////////////////////////////////
    // RENDER SCENE
    /////////////////////////////////////////////////////////////////////

    // update shadow maps (if any)
    world->updateShadowMaps(false, mirroredDisplay);

    // render world
    camera->renderView(width, height);

    // wait until all OpenGL commands are completed
    glFinish();

    // check for any OpenGL errors
    GLenum err;
    err = glGetError();
    if (err != GL_NO_ERROR) cout << "Error:  %s\n" << gluErrorString(err);
}

//------------------------------------------------------------------------------

void updateHaptics(void)
{
    // simulation in now running
    simulationRunning  = true;
    simulationFinished = false;

    cPrecisionClock clock;
    clock.start(true);

    std::string filename = "/home/aldo/DiffInvertedPendulum/trial.csv";

    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file '" << filename << "' for writing." << std::endl;
        return;
    }

    double simTime = 0;

    // main haptic simulation loop
    while(simulationRunning)
    {

        // time the loop
        clock.stop();
        double dt = clock.getCurrentTimeSeconds();
        simTime += dt;
        clock.start(true);

        device->updateFromDevice();

        // update the dynamics of the simulation
        pendulum->ImplicitEuler(dt);

        // send computed force, torque, and gripper force to haptic device
        hapticDevice->setForceAndTorqueAndGripperForce(-cVector3d(pendulum->Fc), -cVector3d(pendulum->Tc), 0);

        writeData(file,simTime);

        // signal frequency counter
        freqCounterHaptics.signal(1);
    }

    // close rhe recording file
    file.close();

    // exit haptics thread
    simulationFinished = true;
}

void writeData(std::ofstream& file, double simTime)
{
    Vector3d x = pendulum->x;
    Vector3d v = pendulum->P/pendulum->m;
    Quaterniond q = pendulum->q;
    Vector3d w = pendulum->q*pendulum->J*pendulum->q.inverse()*pendulum->L;

    Vector3d xh = pendulum->xh;
    Vector3d vh = pendulum->vh;
    Quaterniond qh = pendulum->qh;
    Vector3d wh = pendulum->wh;

    Vector3d f = pendulum->Fc;
    Vector3d t = pendulum->Tc;

    // Writing the data to the file
    file << simTime << ","
         << x(0) << "," << x(1) << "," << x(2) << ","
         << v(0) << "," << v(1) << "," << v(2) << ","
         << q.w() << "," << q.x() << "," << q.y() << "," << q.z() << ","
         << w(0) << "," << w(1) << "," << w(2) << ","
         << xh(0) << "," << xh(1) << "," << xh(2) << ","
         << vh(0) << "," << vh(1) << "," << vh(2) << ","
         << qh.w() << "," << qh.x() << "," << qh.y() << "," << qh.z() << ","
         << wh(0) << "," << wh(1) << "," << wh(2) << ","
         << f(0) << "," << f(1) << "," << f(2) << ","
         << t(0) << "," << t(1) << "," << t(2)
         << std::endl;
}
