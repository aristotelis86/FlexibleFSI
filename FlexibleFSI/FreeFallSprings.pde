/**********************************************************************
      FreeFallSprings class: Creates a demo to simulate
      multiple springs in free-fall. Point-Point collision
      is active.

Example code:
int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution
int N = 6;
PVector gravity = new PVector(0,3);
Window view; // convert pixels to non-dim frame
FreeFallSprings demo;

void settings(){
    size(600, 600);
}

void setup() {  
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  demo = new FreeFallSprings( N, gravity, view );
} // end of setup

void draw() {
  background(185);  
  demo.RunDemo();
}
**********************************************************************/

class FreeFallSprings {
  
  //================= Attributes ====================//
  int N;
  PVector gravity;
  Window view;
  CollisionHandler collider;
  ControlPoint [] cpoints;
  Spring [] springs;
  
  //================= Constructor ====================//
  FreeFallSprings( int N_, PVector f, Window view_ ) {
    N =N_;
    gravity = f;
    view = view_;
    
    cpoints = new ControlPoint[2*N];
    springs = new Spring[N];
    
    for (int i=0; i<N; i++) {
      cpoints[2*i] = new ControlPoint( new PVector(view.x.inE/4.,view.y.inE*(i+1)/(float)N), 5,  10, view );
      cpoints[2*i+1] = new ControlPoint( new PVector(3*view.x.inE/4.,view.y.inE*(i+1)/(float)N), 5,  10, view );
      float th = (cpoints[2*i].diameter + cpoints[2*i+1].diameter)/2.;
      springs[i] = new Spring(cpoints[2*i], cpoints[2*i+1], (i+1)*nx/8., 5, 0.5, th, view );
    }
    
    collider = new CollisionHandler( cpoints, springs );
  }
  
  //================= Methods =====================//
  void RunDemo() {
    for (ControlPoint cp : cpoints) cp.clearForce();
    for (Spring sp: springs) sp.ApplyAllForces();
    for (ControlPoint cp : cpoints) {
      cp.ApplyForce( gravity );
      cp.updateAlt( 0.1 );
      cp.updateAlt2( 0.1 );
    }
    collider.HandleCollisions();
    for (ControlPoint cp : cpoints) cp.display();
    for (Spring sp: springs) sp.display();
  }
  
}