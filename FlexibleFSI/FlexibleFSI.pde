/************************ Input Section ************************/

int nx = 50; // x-dir resolution
int ny = 50; // y-dir resolution

float x, y, vx, vy;
/************************ Setup Section ************************/
Window view; // convert pixels to non-dim frame
ControlPoint cpoints;

void settings(){
    size(450, 450);
}

void setup() {
  
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  cpoints = new ControlPoint( new PVector(nx/2., ny/2. ), 10,  5, view );
  x = nx/2.;
  y = ny/2.;
  vx = 0; 
  vy = 0;
  println(view.x.outS);
  println(view.x.outE);
} // end of setup


/************************ Draw Section ************************/
void draw() {
  background(185);
  
  cpoints.display();
  cpoints.randomWalk();
}