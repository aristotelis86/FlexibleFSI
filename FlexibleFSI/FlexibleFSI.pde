/************************ Input Section ************************/

int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution
int N = 15;


PVector gravity = new PVector(6,6);

/************************ Setup Section ************************/
Window view; // convert pixels to non-dim frame
SimpleCollisionTest Demo;

void settings(){
    size(600, 600);
}

void setup() {
  
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  
  Demo = new SimpleCollisionTest( N, gravity, view );
  
} // end of setup


/************************ Draw Section ************************/
void draw() {
  background(185);
  
  Demo.runDemo();
}


void keyPressed() {
  Demo.terminateDemo();
}