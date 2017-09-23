/************************ Input Section ************************/

int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution

float x, y, vx, vy;
PVector gravity = new PVector(5,-5);

PrintWriter myInfo;
/************************ Setup Section ************************/
Window view; // convert pixels to non-dim frame
ControlPoint cpoints;

void settings(){
    size(450, 450);
}

void setup() {
  
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  cpoints = new ControlPoint( new PVector(nx/2., ny/2. ), 5,  10, view );
  //cpoints.makeFixedx();
  
  myInfo = createWriter("./info/cpoints.txt");
} // end of setup


/************************ Draw Section ************************/
void draw() {
  background(185);
  
  cpoints.display();
  cpoints.clearForce();
  cpoints.ApplyForce( gravity );
  cpoints.updateAlt( 0.1 );
  cpoints.updateAlt2( 0.1 );
  cpoints.BoundCollision(.5);
  cpoints.dampInfo(myInfo);
}


void keyPressed() {
  myInfo.flush();  // Writes the remaining data to the file
  myInfo.close();  // Finishes the file
  exit();  // Stops the program
}