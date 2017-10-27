//int nx = 150; // x-dir resolution
//int ny = 150; // y-dir resolution
//int N = 15;
//PVector gravity = new PVector(0,1);
//Window view; // convert pixels to non-dim frame
//SimpleCollisionTest Demo;

//void settings(){
//    size(600, 600);
//}

//void setup() {
//  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
//  Demo = new SimpleCollisionTest( N, gravity, view );
//} // end of setup

//void draw() {
//  background(185); 
//  Demo.runDemo();
//}

//void keyPressed() {
//  Demo.terminateDemo();
//}

int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution
int N = 3; // number of control points to create 

Window view; // convert pixels to non-dim frame
ControlPoint [] cpoints = new ControlPoint[N];

void settings(){
    size(600, 600);
}

void setup() {
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  for (int i=0; i<N; i++) {
    cpoints[i] = new ControlPoint( new PVector(random(view.x.inE),random(view.y.inE)), 5,  5, view );
  }
}

void draw() {
  for (ControlPoint cp : cpoints) {
    cp.randomWalk();
    cp.display();
  }
}