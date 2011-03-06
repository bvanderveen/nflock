import fullscreen.*;

class NFlock
{
    PApplet pApplet;
    
    class Bird
    {
        public double[] p, v, p0;
        public Bird(double[] p, double[] v) { this.p = p; this.v = v; }
    }
    
    int 
        count = 200,
        dimension = 6,
        world_width = 1280,
        world_height = 720,
        world_depth = 1280;

    double
        phase = 0,
        phaseRate = .002,
        visibility = .000001,
        cohesion = .05,
        alignment = .005,
        separation = .001,
        per_bird_separation = .003,
        attractor_strength = .0003,
        momentum = .002,
        randomness = .001,
        stall = 1;

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
        attractor = new double[] { .5, .5, .5, .5, .5, .5 },
        attractor1 = new double[] { .5, .5, .5, .5, .5, .5 },
        random_min = new double[] { -.02, -.02, -.02, -.01, -.01, -.01},
        random_max = new double[] { .02, .02, .02, .01, .01, .01 };
        

    Bird[] birds;
    
    NFlock initialize()
    {
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
        int num_birds_in_sight = 0;
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

            num_birds_in_sight++;

            // arithmetic mean
            for (int j = 0; j < dimension; j++)
                //cohesion_point[j] = (num_birds_in_sight * cohesion_point[j] + birds[i].p[j]) / (num_birds_in_sight + 1);
                cohesion_point[j] += birds[i].p[j];
            
            // simple sum, normalize later
            for (int j = 0; j < dimension; j++)
                alignment_vector[j] += birds[i].v[j];
                
            //println("alignment");
            //printArray(alignment_vector);

            // weighted sum (weight proportional to square of proximity)
            double d2 = distance_squared(bird.p, birds[i].p);
            //println("distance squared: " + d2);
            for (int j = 0; j < dimension; j++)
                separation_vector[j] += (per_bird_separation * (bird.p[j] - birds[i].p[j])) / d2;
                //separation_vector[j] += bird.p[j] - birds[i].p[j];
        }
        
        
        //attractor[0] = mouseX;
        //attractor[1] = mouseY;
        attractor[0] = (double)mouseX / (double)world_width;
        attractor[1] = (double)mouseY / (double)world_height;
        //attractor[2] = ((double)(mouseX + mouseY) / (world_width + world_height));
        attractor[3] = .5 * Math.sin(phase) + .5;
        attractor[4] = .5 * Math.sin(phase * phase * .8 + 1) + .5;
        attractor[5] = .5 * Math.sin(phase + 2) + .5;
        
        attractor1[0] = (double)(world_width - mouseX) / (double)world_width;
        attractor1[1] = (double)(world_height - mouseY) / (double)world_height;
        //attractor[2] = ((double)(mouseX + mouseY) / (world_width + world_height));
        attractor1[3] = .5 * Math.sin(phase * phase * .8) + .5;
        attractor1[4] = .5 * Math.sin(phase + 1) + .5;
        attractor1[5] = .5 * Math.sin(phase - 1) + .5;
        
        //println("solved with " + num_birds_in_sight + " birds in sight");
        double[] to_attractor = difference(attractor, bird.p);
        double[] to_attractor1 = difference(attractor1, bird.p);
        //double attractor_distance_squared = distance_squared(to_attractor, attractor);
        //double[] to_cp = difference(cohesion_point, bird.p);
        //double cp_distance_squared = distance_squared(cohesion_point, to_cp);
        
        bird.v = sum(new double[][] {
            scale(cohesion, normalize(difference(cohesion_point, bird.p))),
            //scale(cohesion * cp_distance_squared, normalize(to_cp)),
            scale(alignment, normalize(alignment_vector)),
            scale(separation, normalize(separation_vector)),
            //scale(1/attractor_strength * Math.sqrt(attractor_distance_squared), normalize(to_attractor)),
            scale(attractor_strength, normalize(to_attractor)),
            scale(attractor_strength, normalize(to_attractor1)),
            scale(momentum, normalize(bird.v)),
            scale(randomness, random(random_min, random_max))
        });
        
        /*
        double[] v = new double[] { bird.v[0], bird.v[1] };
        double m2 = magnitude_squared(v);
        
        double r = m2 / stall;
        if(r < 1)
        {
            double[] rv = scale(1 / Math.sqrt(m2), v);
            bird.v[0] = rv[0];
            bird.v[1] = rv[1];
        }*/
        
        //double m2 = magnitude_squared(bird.v);
        
        //double r = m2 / stall;
        //if(r < 1)
        //    bird.v = scale(1 / Math.sqrt(m2), bird.v);
    }

    void draw_bird(Bird bird)
    {
        bird.p0 = bird.p;      
        bird.p = sum(new double[][] { bird.p, bird.v });
       
        /*
        tline(
          (int)bird.p0[0],
          (int)bird.p0[1],
          (int)bird.p[0],
          (int)bird.p[1],
          (int)bird.p[3],
          (int)bird.p[4],
          (int)bird.p[5],
          (float)(bird.p[2] / world_depth) * .1);
        */
        /*
        color c = color((int)bird.p[3],
          (int)bird.p[4],
          (int)bird.p[5],
          (float)(bird.p[2] / world_depth) * 255);
        stroke(c);
        line(
          (int)bird.p0[0] + (bird.p0[0] < bird.p[0] ? 1 : -1),
          (int)bird.p0[1] + (bird.p0[1] < bird.p[1] ? 1 : -1),
          (int)bird.p[0],
          (int)bird.p[1]);
          */
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
int iterations = 0;
void setup() { 
  noCursor();
  flock = new NFlock().initialize(); 
  fs = new FullScreen(this);
  fs.enter();
}
void draw() { if(keyPressed || true) for (int i = 0; i < 5; i++) flock.iterate(); }

/*void setup()
{
  double[] i = new double[1];
  println("it is" + i[0]);
}*/
