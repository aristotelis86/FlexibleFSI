/************************ Input Section ************************/

int nx = 50; // x-dir resolution
int ny = 50; // y-dir resolution

/************************ Setup Section ************************/
Window view; // convert pixels to non-dim frame
ControlPoint cpoints;

void settings(){
    size(450, 450);
}

void setup() {
  
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  cpoints = new ControlPoint( new PVector(nx/2., ny/2. ), 10,  49, view );
  
  
} // end of setup


/************************ Draw Section ************************/
void draw() {
  cpoints.display();
  //noLoop();
}