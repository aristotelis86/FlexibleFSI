int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution

float Length = 20;
float thick = 1;
float MassNum = 0.1;
int resol = 3;
float stiffness = 100;
PVector lpos = new PVector(nx/4.,ny/3.);
PVector align = new PVector(4,0);
FlexibleSheet sheet;
Window view; // convert pixels to non-dim frame

PVector gravity = new PVector(0,5);
float t = 0;
float dt;
float maxVel = 20;

void settings(){
    size(600, 600);
}

void setup() {
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  sheet = new FlexibleSheet( Length, thick, MassNum, resol, stiffness, lpos, align, view );
  sheet.cpoints[0].makeFixed();
  
  sheet.Calculate_Stretched_Positions( gravity );
  // Create the distortion
  int N = sheet.numOfpoints;
  
  // Add an impulse (x-dir) to the particles
  for (int i = 1; i < N; i++) {
    sheet.cpoints[i].velocity.x += ((i-1)/(N-2)) * maxVel;
  }
  dt = sheet.dtmax;
} // end of setup

void draw() {
  background(185); 
  
  sheet.updateAlt( dt, gravity );
  sheet.updateAlt2( dt, gravity );
  
  sheet.display(color(0, 255, 0));
}