/************************ Input Section ************************/

int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution
int N = 2;

float x, y, vx, vy;
PVector gravity = new PVector(25,25);

/************************ Setup Section ************************/
Window view; // convert pixels to non-dim frame
ControlPoint [] cpoints;
PrintWriter [] myInfo;
CollisionHandler collider;

void settings(){
    size(600, 600);
}

void setup() {
  
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  cpoints = new ControlPoint[N];
  myInfo = new PrintWriter[N];
  
  for (int i=0; i<N; i++) {
    cpoints[i] = new ControlPoint( new PVector(random(nx),random(ny)), 5,  10, view );
    myInfo[i] = createWriter("./info/cpoints"+i+".txt");
  }
  
  cpoints[0].UpdatePosition(nx,ny);
  cpoints[1].UpdatePosition(nx/2.,ny/2.);
  collider = new CollisionHandler( cpoints );
  
} // end of setup


/************************ Draw Section ************************/
void draw() {
  background(185);
  
  for (ControlPoint cp : cpoints) {
    cp.clearForce();
    cp.ApplyForce( gravity );
    cp.updateAlt( 0.1 );
    cp.updateAlt2( 0.1 );
  }
  
  collider.HandleCollisions();
  for (ControlPoint cp : cpoints) cp.display();
  for (int i=0; i<N; i++) cpoints[i].dampInfo(myInfo[i]);
}


void keyPressed() {
  for (int i=0; i<N; i++) {
    myInfo[i].flush();  // Writes the remaining data to the file
    myInfo[i].close();  // Finishes the file
  }
  exit();  // Stops the program
}