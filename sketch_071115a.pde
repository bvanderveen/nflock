import fullscreen.*;

class NFlock
{
    PApplet pApplet;
    
    class Bird
    {
        public double[] p, v, p0;
        public Bird(double[] p, double[] v) { this.p = p; this.v = v; }
    }
    
    public int 
        count = 200,
        dimension = 6,
        world_width,
        world_height,
        world_depth = 1280;

    double
        phase = 0,
        phaseRate = .002,
        visibility = .0005,
        cohesion = .000007,
        alignment = .0005,
        separation = 0.001,
        per_bird_separation = .003,
        drive = .003,
        current_drive = 0;
        first_attractor_strength = .00000,
        first_attractor_strength_pressed = .004,
        attractor_strength = .001,
        momentum = .002,
        randomness = .000,
        stall = 1;

    public double[][] attractors;

    double[] 
        //init_pmax = new double[] { world_width - (world_width * .8), world_height - (world_height * .9), world_depth - (world_depth * .9), 255, 255, 255 },
        init_pmax = new double[] { 1, 1, 1, 1, 1, 1 },
        //init_pmin = new double[] { world_width * .8, world_height * .8, world_depth * .0, 0, 0, 0 },
        init_pmin = new double[] { 0, 0, 0, 0, 0, 0 },
        //init_vmax = new double[] { -1, .3, -1, 0.0015, 0.0025, 0.0035 }, 
        init_vmax = new double[] { .05, .05, .05, 0.1, 0.1, 0.1 }, 
        //init_vmin = new double[] { -.8, .4, -2, -.001, -.0090, -.0070 },
        init_vmin = new double[] { -.05, -.05, -.05, -.1, -.1, -.1 },
        //attractor = new double[] { world_width / 2, world_height / 2, world_depth - (world_depth *  .6), 20, 118, 107 },
        random_min = new double[] { -.02, -.02, -.02, -.01, -.01, -.01},
        random_max = new double[] { .02, .02, .02, .01, .01, .01 };
        

    Bird[] birds;
    
    NFlock initialize(int w, int h)
    {
        world_width = w;
        world_height = h;
        size(world_width, world_height);
        background(#000000);
        
        // fix values
        visibility = Math.pow(visibility,2);
        stall = Math.pow(stall,2);
        
        birds = new Bird[count];

        for (int i = 0; i < birds.length; i++)
            birds[i] = new Bird(random(init_pmax, init_pmin),  random(init_vmax, init_vmin));
         
         return this;
    }
    
    void iterate()
    {
      phase += phaseRate;
      if (phase > Math.PI * 2)
        phase = 0;
      for (int j = 0; j < birds.length; j++)
      {
          solve_bird(birds[j],j);
          draw_bird(birds[j]);
      }
      
      //trect(0,0,world_width,world_height,255,255,255,.5);
      fill(color(0,0,0,1));
      rect(0,0,world_width,world_height);
      //background(#ffffff);
    }
    
    void solve_bird(Bird bird, int index)
    {
        //println("solving");
        double[]
            cohesion_point = new double[dimension],
            alignment_vector = new double[dimension],
            separation_vector = new double[dimension];

        System.arraycopy(bird.p, 0, cohesion_point, 0, dimension);
        System.arraycopy(bird.v, 0, alignment_vector, 0, dimension);

        for (int i = 0; i < count; i++)
        {
            if(i == index)
                continue;
                
            if (distance_squared(bird.p, birds[i].p) > visibility)
                continue;

            for (int j = 0; j < dimension; j++)
                cohesion_point[j] += birds[i].p[j];
            
            // simple sum, normalize later
            for (int j = 0; j < dimension; j++)
                alignment_vector[j] += birds[i].v[j];
            
            for (int j = 0; j < dimension; j++)
                separation_vector[j] += bird.p[j] - birds[i].p[j];
        }
        
        double[] driver = new double[dimension];
        
        driver[0] = (double)mouseX / (double)world_width;
        driver[1] = (double)mouseY / (double)world_height;
        //attractor[2] = ((double)(mouseX + mouseY) / (world_width + world_height));
        driver[2] = .5;
        driver[3] = .5 * Math.sin(phase) + .5;
        driver[4] = .5 * Math.sin(phase * phase * .8 + 1) + .5;
        driver[5] = .5 * Math.sin(phase + 2) + .5;
        
        double[][] components = new double[6 + attractors.length][];
        
        components[0] = scale(cohesion, normalize(difference(cohesion_point, bird.p)));
        components[1] = scale(alignment, normalize(alignment_vector));
        components[2] = scale(separation, normalize(separation_vector));
        components[3] = scale(momentum, normalize(bird.v));
        components[4] = scale(randomness, random(random_min, random_max));
        components[5] = scale(drive, normalize(difference(driver, bird.p)));
        double as = attractor_strength / attractors.length;
        for (int i = 0; i < attractors.length; i++)
        {
          components[i + 6] = scale(as, normalize(difference(attractors[i], bird.p)));
        }
        
        bird.v = sum(components);
    }

    void draw_bird(Bird bird)
    {
        bird.p0 = bird.p;      
        bird.p = sum(new double[][] { bird.p, bird.v });
       
        color c = color((int)bird.p[3] * 255,
          (int)(bird.p[4] * 255),
          (int)(bird.p[5] * 255),
          (int)(bird.p[2] * 255));
        stroke(c);
        
        int p0x = (int)(bird.p0[0] * world_width);
        int px = (int)(bird.p[0] * world_width);
        int p0y = (int)(bird.p0[1] * world_height);
        int py = (int)(bird.p[1] * world_height);
        
        line(
          (int)p0x,// + (p0x < px ? 1 : -1),
          (int)p0y,// + (p0y < py ? 1 : -1),
          (int)px,
          (int)py);
          
        //println("drawing bird: " + px + " " + py);
        
    }

    double[] random(double[] upper, double[] lower)
    {
        double[] r = new double[upper.length];
        for (int i = 0; i < r.length; i++)
            r[i] = lower[i] + Math.random() * (upper[i] - lower[i]);
        return r;
    }

    double distance_squared(double[] a, double[] b)
    {
        double r = 0, t = 0;
        for (int i = 0; i < a.length; i++)
        {
            t = (a[i] - b[i]);
            r += t * t;
        }
        return r;
    }

    double magnitude_squared(double[] a)
    {
        return distance_squared(a, new double[a.length]);
    }

    double[] sum(double[][] a)
    {
        double[] r = new double[a[0].length];
        for (int i = 0; i < a.length; i++)
            for (int j = 0; j < a[0].length; j++)
                r[j] += a[i][j];
        return r;
    }

    double[] difference(double[] a, double[] b)
    {
        double[] r = new double[a.length];
        for (int i = 0; i < a.length; i++)
            r[i] = a[i] - b[i];
        return r;
    }

    double[] scale(double a, double[] x)
    {
        double[] r = new double[x.length];
        for (int i = 0; i < x.length; i++)
            r[i] = a * x[i];
        return r;
    }

    double[] normalize(double[] a)
    {
        //println("normalizing");
        //printArray(a);
        // slow. better way to find invsqrt?
        
        double m = magnitude_squared(a);
        
        if(Math.abs(m) < .0000000000001)
          m = .0000000000001;
        return scale(1 / Math.sqrt(m), a);
    }   
}

NFlock flock;
FullScreen fs;
ArrayList<double[]> attractorSetup;
boolean done;

void setup() { 
  flock = new NFlock().initialize(800,800); 
  //fs = new FullScreen(this);
  //fs.enter();
  attractorSetup = new ArrayList<double[]>();
}

void mousePressed()
{
  if (done)
    return;
    
  double[] a = new double[flock.dimension];
  a[0] = ((double)mouseX / flock.world_width);
  a[1] = ((double)mouseY / flock.world_height);
  a[2] = .5;
  a[3] = .5;
  a[4] = .5;
  a[5] = .5;
  
  attractorSetup.add(a);
  fill(255);
  ellipse(mouseX - 2, mouseY - 2, 4, 4);
}

void keyReleased()
{
  if (key == 'f')
    flock.first_attractor_strength = .00000000000001;
}

void keyPressed()
{
  if (key == 'f')
    flock.first_attractor_strength = flock.first_attractor_strength_pressed;
   
  if (key == 'd')
  {
    flock.attractors = new double[attractorSetup.size()][];
    
    for (int i = 0; i < flock.attractors.length; i++)
    {
      flock.attractors[i] = attractorSetup.get(i);
    }
    done = true;
    print("done");
    noCursor();
  }
}

void draw() { 
  if (done)
  {
    for (int i = 0; i < 5; i++) flock.iterate(); 
  }
}

/*void setup()
{
  double[] i = new double[1];
  println("it is" + i[0]);
}*/
