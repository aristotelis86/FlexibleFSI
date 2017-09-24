/************************ Input Section ************************/

int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution
int N = 2;


PVector gravity = new PVector(0,0);

/************************ Setup Section ************************/
Window view; // convert pixels to non-dim frame
ControlPoint [] cpoints;
Spring spring;
CollisionHandler collider;

void settings(){
    size(600, 600);
}

void setup() {
  
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  
  cpoints = new ControlPoint[N];
  cpoints[0] = new ControlPoint( new PVector(nx/3.,ny/2.), 5,  10, view );
  cpoints[1] = new ControlPoint( new PVector(2*nx/3.,ny/2.), 5,  10, view );
  
  spring = new Spring(cpoints[0], cpoints[1], nx/8., 5, 0.5, 1, view );
  
  collider = new CollisionHandler( cpoints );
} // end of setup


/************************ Draw Section ************************/
void draw() {
  background(185);
  for (ControlPoint cp : cpoints) cp.clearForce();
  spring.ApplyAllForces();
  for (ControlPoint cp : cpoints) {
    cp.ApplyForce( gravity );
    cp.updateAlt( 0.1 );
    cp.updateAlt2( 0.1 );
  }
  collider.HandleCollisions();
  for (ControlPoint cp : cpoints) cp.display();
  spring.display();
  //saveFrame("./movie/frame_######.png");
}


void keyPressed() {
  //Demo.terminateDemo();
}