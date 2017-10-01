int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution

float Length = 20;
float thick = 1;
float MassNum = 0.1;
int resol = 3;
float stiffness = 100;
PVector lpos = new PVector(nx/4.,ny/3.);
PVector align = new PVector(4,0);
FlexibleSheet [] sheet;
Window view; // convert pixels to non-dim frame
WriteInfo writer;

PVector gravity = new PVector(0,5);
float t = 0;
float dt;
float maxVel = 20;

void settings(){
    size(600, 600);
}

void setup() {
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  sheet = new FlexibleSheet[2];
  sheet[0] = new FlexibleSheet( Length, thick, MassNum, resol, stiffness, lpos, align, view );
  sheet[0].cpoints[0].makeFixed();
  
  sheet[0].Calculate_Stretched_Positions( gravity );
  // Create the distortion
  int N = sheet[0].numOfpoints;
  
  // Add an impulse (x-dir) to the particles
  for (int i = 1; i < N; i++) {
    sheet[0].cpoints[i].velocity.x += ((i-1)/(N-2)) * maxVel;
  }
  sheet[1] = new FlexibleSheet( Length, thick, 1.5*MassNum, resol, 3*stiffness, lpos.add(new PVector(30,0)), align, view );
  sheet[1].cpoints[0].makeFixed();
  
  sheet[1].Calculate_Stretched_Positions( gravity );
  // Create the distortion
  N = sheet[1].numOfpoints;
  
  // Add an impulse (x-dir) to the particles
  for (int i = 1; i < N; i++) {
    sheet[1].cpoints[i].velocity.x += ((i-1)/(N-2)) * maxVel;
  }
  
  dt = sheet[0].dtmax;
  if (sheet[1].dtmax<dt) dt=sheet[1].dtmax;
  
  writer = new WriteInfo(sheet);
} // end of setup

void draw() {
  background(185); 
  
  for (FlexibleSheet fs : sheet) {
    fs.updateAlt( dt, gravity );
    fs.updateAlt2( dt, gravity );
    fs.display(color(0, 255, 0));
  }
  
  writer.dampAllInfo( t, true );
  
  t += dt;
}


void keyPressed() {
  writer.terminateFiles();
}