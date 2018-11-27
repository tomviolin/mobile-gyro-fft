package processing.test.gyroscope;

import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import android.os.Bundle; 
import ketai.sensors.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class Gyroscope extends PApplet {





KetaiSensor sensor;
float rotationX, rotationY, rotationZ;
int NR=512;
float pi = atan(1)*4;
float kdt=1/420.f;
float kf = .4f;
float v[];
float d[];
float e[];
float ep[], es[], ke[];
float savnum=0;
int om=0;
public float x2f(float x) {
  return pow(1.005f, x)/pow(1.005f, NR)*18-0.3f;
}
public float x2m(float x) {
  return 3;//20*exp(-x/NR);///x2f(x); //*Math.exp(-x/1000);//0./(x+10); //exp(-x)/exp(-256)*10.;
}
public float x2d(float x) {
  return 1;
}
public void mousePressed() {
  // if (/ != null) mic.start();
  om=(om+1)%4;
}

class CircBuff {

  float[] buff;
  int bufrd=0;
  int bufwr=0;
  int size=0;

  CircBuff(int siz) {
    size=siz;
    buff = new float[size];
  }

  public void write(float x) {
    int nextWrite = (bufwr+1)%size;
    if (nextWrite == bufrd) {
      // buffer full, throw away oldest value
      read();
    }
    buff[nextWrite] = x;
    bufwr=nextWrite;
  }
  public float read() {
    if (bufwr == bufrd) {
      return Float.NaN;
    } else {
      float val=buff[bufrd];
      bufrd = (bufrd+1)%size;
      return val;
    }
  }
  public int left() {
    // how many left to read?
    if (bufrd == bufwr) return 0;
    if (bufrd < bufwr) return bufwr-bufrd;
    return bufwr+size-bufrd;
  }
  public float[] readall() {
    if (left()==0) return null;
    float[] ans = new float[left()];
    if (bufrd <  bufwr) {
      arrayCopy(buff, bufrd, ans, 0, left());
      bufrd=bufwr;
    } else {
      arrayCopy(buff, bufrd, ans, 0, size-bufrd);
      arrayCopy(buff, 0, ans, size-bufrd, left()-(size-bufrd));
      bufrd=bufwr;
    }
    return ans;
  }
}
int wfy;
PShader scroll, scrolltop;
boolean oktodraw=false;
CircBuff xbuff;

public void setup()
{
  
  background(0);
  frameRate(60);
  rectMode(CORNER);
  xbuff = new CircBuff(10000);
  sensor = new KetaiSensor(this);
  sensor.setSamplingRate(0);

  orientation(PORTRAIT);
  //android.app.Activity.getWindow().addFlags(WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);

  println("continuing...");
  wfy=height/4;
  d = new float[NR];
  e = new float[NR];
  v = new float[NR];
  ke= new float[NR];
  ep= new float[NR];
  es= new float[NR];
  for (int i = 0; i<NR; ++i) {
    v[i]=0.0f;
    ke[i]=Float.NaN;
  }
  textAlign(CENTER, CENTER);
  PrintWriter output = createWriter("scroll.glsl");
  output.println("\n"+
    "precision highp float;\n"+
    "#define PROCESSING_TEXTURE_SHADER\n"+
    "uniform vec2 texSize;\n"+
    "uniform sampler2D texture;\n"+
    "void main(void){\n"+
    "vec2 TM=1./texSize;\n"+
    "vec2 UV=gl_FragCoord.xy*TM;\n"+
    "if (UV.y > 3./4.) { \n" +
    "   gl_FragColor = texture2D(texture,UV);\n"+
    "}\n"+
    " else  { "+
    "   gl_FragColor = texture2D(texture, UV+TM*vec2(0.,1.0));\n"+
    " }}");
  output.flush();
  output.close();
  scroll=loadShader("scroll.glsl");
  scroll.set("texSize", PApplet.parseFloat(width), PApplet.parseFloat(height));

  output = createWriter("scrolltop.glsl");
  output.println("\n"+
    "precision highp float;\n"+
    "#define PROCESSING_TEXTURE_SHADER\n"+
    "uniform vec2 texSize;\n"+
    "uniform sampler2D texture;\n"+
    "uniform float h;\n"+
    "void main(void){\n"+
    "vec2 TM=1./texSize;\n"+
    "vec2 UV=gl_FragCoord.xy*TM;\n"+
    "if (UV.y > 7./8. && UV.x > TM.x*h-1.5) { \n" +
    "   gl_FragColor = texture2D(texture,UV+TM*vec2(h,0.));\n"+
    "}\n"+
    " else  { "+
    "   gl_FragColor = texture2D(texture, UV);\n"+
    " }}");
  output.flush();
  output.close();
  scrolltop=loadShader("scrolltop.glsl");
  scrolltop.set("texSize", PApplet.parseFloat(width), PApplet.parseFloat(height));
  scrolltop.set("h", 1.0f);
  println(width+"x"+height);
  oktodraw=true;
  sensor.start();
}

float noisefloor=0;
// this trick keeps the screen on.
public void onCreate(Bundle bundle) {
  println("creating...");
  super.onCreate(bundle);
  getActivity().getWindow().addFlags(android.view.WindowManager.LayoutParams.FLAG_KEEP_SCREEN_ON);
}
int graphX=0;
public void draw()
{
  if (kf>.0005f) {
    background(0);
    fill(255, 150, 0);
    textSize(height/32);
    text("CALIBRATING... "+nf(kf, 1, 5), width/2, height/2);
    if (kf < .001f) {
      for (int i=0; i<NR; ++i) {
        if (e[i] > noisefloor) noisefloor=e[i];
      }
    }
    return;
  }
  if (!oktodraw) return;
  float nextfreq=1.0f;// next freq tick mark
  float nextftik=1.0f; // minor ticks
  // reset freq response graph field
  fill(0, 120, 255);
  stroke(0);
  rect(0, height/8, width, height/8);
  //  rect(0,0,width,height/8);
  /*
  int SPACING=10;
   float[] xr = xbuff.readall();
   if (xr != null){
   println("got "+xr.length+" readings!");
   for (int i=0; i< xr.length-1; ++i){
   float graphy1 = atan((float)(xr[i  ]*50.))/6.28*height/8.;
   float graphy2 = atan((float)(xr[i+1]*50.))/6.28*height/8.;
   stroke(color(255, 200, 0));
   strokeWeight(5);
   line((float)(width-(xr.length-i*SPACING)), height/16+(float)graphy2, (float)(width-(xr.length-(i-1)*SPACING)), height/16+(float)graphy1);
   }
   }
   */
  //for (int i=0;i<NR;++i){
  //  rect((float)(width*i/NR),(float)(height/16+d[i]),(float)(width/NR),1);
  //}
  //filter(scrolltop);
  stroke(255);
  fill(255);

  float maxe = e[1];
  float mine = e[1];
  float tote = 0;
  for (int i = NR/32; i<NR-1; ++i) {
    if (Float.isNaN(ke[i])) ke[i]=e[i];
    else {
      ke[i] = ke[i]*.9f+e[i]*.1f;
      if (maxe < ke[i]) maxe = ke[i];
      if (mine > ke[i]) mine = ke[i];
      tote += ke[i];
      if (i>0) ep[i]=ke[i]-ke[i-1];
      if (i>1) es[i]=ep[i]-ep[i-1];
    }
  }
  //=maxe;
  tote /= NR;

  float maxdisp = tote*3;
  if (maxdisp < maxe) maxdisp=maxe;
  if (maxdisp < noisefloor*2) maxdisp=noisefloor*2;
  float boo = ((tote-mine)/(maxe-mine));



  for (int i = NR/32; i<NR-1; ++i) {
    //float m=x2m(i);
    float f=x2f(i);
    //float k = m*(f*2*pi)*(f*2*pi);
    boolean ticked=false;
    float er=(ke[i])/maxdisp;
    if (f>=nextfreq) {
      fill(255);
      stroke(0);
      nextfreq*=2;
      strokeWeight(5);
      line(i*width/NR, height/4, 
        i*width/NR, height/8+25);
      textSize(25);
      text(nf(((float)nextfreq/2), 1, 0), (float)(i*width/NR), (float)(height/7.5f));
      ticked=true;
    } else {
      fill(255);
      stroke(255);
      rect((float)(i*width/NR), (float)(height/4), (float)(ceil((float)width/(float)NR)), (float)(-min(1.f, (float)er)*height/16.f));//+d[i]*0*10000.);
    }
    if (f>=nextftik) {
      fill(0);
      stroke(0);
      strokeWeight(3);
      nextftik++;
      if (!ticked) {
        line(i*width/NR, height/8, 
          i*width/NR, height/8+20);
        textSize(15);
        //textFont(createFont("Arial Black",20));
        text (nf(nextftik-1, 1, 0), i*width/NR, height/8+35);
      }
    }

    if (i>3 && ep[i-1] > 0 && ep[i] <= 0 &&
      ke[i] > maxdisp*.8f) {
      stroke(255, 128, 0);
      strokeWeight(3);
      line(i*width/NR, height*3/16+10, i*width/NR, height*3/16-10);
      textSize(25);
      text(nf(x2f(i), 1, 2), i*width/NR, 
        height*3/16-20);
    }

    strokeWeight(1);

    switch(om) {
    case 0:
      stroke(color((float)er*255, (float)(er*er*255), (float)Math.sqrt(er)*255));
      break;
    case 1:
      stroke((float)er*255);
      break;
    case 2:
      stroke((float)(er*er*255), (float)er*255, (float)Math.sqrt(er)*255);
      break;
    case 3:
      stroke(er*er*255, sqrt(er)*255, er*255);
      break;
    }

    strokeWeight(1);
    rect(i*width/NR, wfy, (i+1)*width/NR-i*width/NR, 0);//+d[i]*0*10000.);
  }
  filter(scroll);
}

float dork=0;
float kx=Float.NaN;
float kkx=Float.NaN;
public void resonate(float xd, float dt) {

  float x;
  if (Float.isNaN(kx)) {
    x=kx=xd;
  } else {
    kx=xd*.5f + kx * .5f;
  }
  if (Float.isNaN(kkx)) {
    kkx=xd;
  } else {
    kkx=xd*.01f+kkx*.99f;
  }

  x= xd-kkx;
  for (int i = 0; i < NR; ++i) {
    //float drag=x2d(i);

    float f = x2f(i);
    float m = 100/f;
    float k = m*(f*2*pi)*(f*2*pi);
    // float k = 10000;
    //float m=k/((f*2*pi)*(f*2*pi));

    //float m = x2m(i);
    float drag = sqrt(m*k)*2/(9);
    //float k=300
    savnum=d[500];
    float v0=v[i];
    float d0=d[i];

    float a0 = -k*d0 - drag*v0 +x/200;      // original acceleration
    float v0_5 = v0 + a0*(dt/2)/m;
    float d0_5 = d0 + v0*(dt/2);
    // compute new full step from halfway
    float a0_5 = -k*d0_5 -drag*v0_5 +x/200;
    float v1_5 = v0_5 + a0_5*(dt)/m;
    float d1_5 = d0_5 + v0_5*(dt);
    // add difference to original position
    float v1 = v0+(a0_5*dt)/m;
    float d1 = d0+(v0_5*dt);
    v[i]=v1;
    d[i]=d1;
    float ep = (d[i]*d[i]*k)/2.f;
    float ek =  m*v[i]*v[i]/2.f;
    e[i]=ek+ep;
  }
}


long lastNano=0;

//void onAccelerometerEvent(float x, float y, float z)
public void onGyroscopeEvent(float x, float y, float z)
{
  //
  //savnum=x;
  rotationX = x;
  rotationY = y;
  rotationZ = z;
  xbuff.write(x);
  long thisNano = System.nanoTime();
  if (lastNano>0) kdt = 1e-9f*(thisNano-lastNano)*kf+kdt*(1-kf);
  // fakery: kdt = 1/420.;
  lastNano=thisNano;
  if (kf>.0005f) {
    kf *=.99f;
    if (kf >=.001f) return;
  }
  resonate(x, kdt);
}
  public void settings() {  fullScreen(P2D); }
}
