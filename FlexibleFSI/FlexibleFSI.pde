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


















//int nx = 150; // x-dir resolution
//int ny = 150; // y-dir resolution
//int N = 2;
//PVector gravity = new PVector(0,5);
//Window view; // convert pixels to non-dim frame
//CollisionHandler collider;
//ControlPoint [] cpoints;
//Spring [] springs;
  
//void settings(){
//    size(600, 600);
//}

//void setup() {  
//  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
//  cpoints = new ControlPoint[2*N];
//  springs = new Spring[N];
  
//  cpoints[0] = new ControlPoint( new PVector(nx/2., ny/10.), 5,  10, view );
//  cpoints[1] = new ControlPoint( new PVector(4*nx/5., 4*ny/10.), 5,  10, view );
//  springs[0] = new Spring(cpoints[0], cpoints[1], 50, 5, 0.1, 1, view );
//  springs[0].matchThickness();
  
//  cpoints[2] = new ControlPoint( new PVector(nx/3., 5*ny/10.), 10,  10, view );
//  cpoints[3] = new ControlPoint( new PVector(2*nx/3., 5*ny/10.), 10,  10, view );
//  springs[1] = new Spring(cpoints[2], cpoints[3], 50, 5, 0.1, 1, view );
//  springs[1].matchThickness();
  
//  collider = new CollisionHandler( cpoints, springs );
//} // end of setup

//void draw() {
//  background(185);  
//  for (ControlPoint cp : cpoints) cp.clearForce();
//  for (Spring sp: springs) sp.ApplyAllForces();
//  for (ControlPoint cp : cpoints) {
//    cp.ApplyForce( gravity );
//    cp.updateAlt( 0.1 );
//    cp.updateAlt2( 0.1 );
//  }
//  collider.HandleCollisions();
//  for (ControlPoint cp : cpoints) cp.display();
//  //for (Spring sp: springs) sp.display();
  
//  //saveFrame("./movie/frame_######.png");
//}