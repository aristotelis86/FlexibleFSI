int nx = 150; // x-dir resolution
int ny = 150; // y-dir resolution

float Length = nx/2.;
float thick = 1;
float MassNum = 1;
int resol = 1;
float stiffness = 10;
PVector lpos = new PVector(nx/4.,ny/3.);
PVector align = new PVector(4,0);
FlexibleSheet sheet;
Window view; // convert pixels to non-dim frame

PVector gravity = new PVector(.1,.1);
void settings(){
    size(600, 600);
}

void setup() {  
  Window view = new Window(1, 1, nx, ny, 0, 0, width, height);
  sheet = new FlexibleSheet( Length, thick, MassNum, resol, stiffness, lpos, align, view );
  
  sheet.Calculate_Stretched_Positions( gravity );
} // end of setup

void draw() {
  background(185); 

  sheet.box.display();
  sheet.display(color(0, 255, 0));
}