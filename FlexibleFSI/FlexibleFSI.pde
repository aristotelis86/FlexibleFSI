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