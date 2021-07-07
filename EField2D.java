import org.opensourcephysics.display.*;
import org.opensourcephysics.display2d.*;
import org.opensourcephysics.displayejs.*;
import javax.swing.*;

public class EField2D {   
    private double q[][];
    private double q_[][];
    private double Nq[][];
    private double nq[][];    
    private double uq[][];
    private double uq_[][];    
    private InteractiveCircle charge[];
    private InteractivePanel panel2D = new InteractivePanel();    
    private double farthestRadius = 0, qRadius = 0.5;
    private int f, g ,h ,i , j, l, m, n, lq, u;
    private int nLines = 14, nLines_, nOrigins, grid = 1;
    private double xis, yp,PSum, x, y, dr, representativeSize = 0.5, originX, originY, rv, R, a, Pxx, Pyy, px, py, Cx, Cy, d, precision, k = 8.987551788 * Math.pow(10, 9), ExSum, EySum;
    private double P[], component[], Dx[], Dy[], Px[], Py[], Ec[], E[], Ex[], Ey[], r[], ExEySum[];
    private javax.swing.JTextField jtext = new javax.swing.JTextField(); 
    private DatasetManager data;   
    private double ang = 0;
    private boolean real = false, trapped = false, allTrapped, trappedLine[];
    
    
    public EField2D(double[][] q_) {        
        q = q_;
        q = new double[q_.length+1][3];
        q[0][0] = 0;
        q[0][1] = 0;
        q[0][2] = 1;
        g = 0;
        while(g < q_.length) {
            q[g+1][0] = q_[g][0];
            q[g+1][1] = q_[g][1];
            q[g+1][2] = q_[g][2];
            g++;
        }
        nq = getNegativeCharges();
        getXCentre(panel2D);
        getYCentre(panel2D);
        eCalc2D(px, py);
        pCalc2D(x,y);
        placeCharges(panel2D);
        //updateChargePosition(xis,yp,panel2D);
        getNegativeCharges();
        checkEntrapment(Pxx,Pyy);
        createEFieldLines2D_N(jtext, precision, nLines, panel2D);
        createEFieldLines2D_P(jtext, precision, nLines, panel2D);
        createEFieldVectors2D(real, grid, representativeSize, panel2D);
        createEFieldIsoPotentialLines2D(panel2D);
    }
    
    public double getXCentre(DrawingPanel panel2D) {
        Cx = panel2D.getXMin() + (panel2D.getWidth() / 2);
        return Cx;
    }
    
    public double getYCentre(DrawingPanel panel2D) {
        Cy = panel2D.getYMin() + (panel2D.getHeight() / 2);
        return Cy;
    }
    
    public double[] eCalc2D(double px, double py) {
        r = new double[q.length];
        E = new double[q.length];
        Ex = new double[q.length];
        Ey = new double[q.length];
        ExEySum = new double[2];
        
        m = 0;
        ExSum = 0;
        EySum = 0;
        while(m < q.length) {
            r[m] = Math.sqrt((px - q[m][0]) * (px - q[m][0]) + (py - q[m][1]) * (py - q[m][1]));
            E[m] = (k * q[m][2]) / (r[m] * r[m]);
            Ex[m] = E[m] * ((px - q[m][0]) / r[m]);
            Ey[m] = E[m] * ((py - q[m][1]) / r[m]);
            ExSum = ExSum + Ex[m];
            EySum = EySum + Ey[m];
            m++;
        }
        
        ExEySum[0] = ExSum;
        ExEySum[1] = EySum;        
        return ExEySum;
    }
    
    public double pCalc2D(double x, double y) {
        r = new double[q.length];
        P = new double[q.length];      
        
        m = 0;
        PSum = 0;
        while(m < q.length) {
            r[m] = Math.sqrt((x - q[m][0]) * (x - q[m][0]) + (y - q[m][1]) * (y - q[m][1]));
            P[m] = (k * q[m][2]) / r[m];            
            PSum = PSum + P[m];            
            m++;
        }
        
        return PSum;
    }
        
    
    public double[][] getNegativeCharges() {
        Nq = new double[q.length-1][3];
        l = 0;
        f = 0;
        while(l < q.length) {
            if(q[l][2] < 0) {
                Nq[f] = q[l];
                f++;                
            }
            
            l++;
        }
        
        return Nq;
    }
    
    public boolean checkEntrapment(double Pxx, double Pyy) {        
        trapped = false;
        lq = 0;        
        while(lq < nq.length) {
            if(Math.sqrt((Pxx - nq[lq][0]) * (Pxx - nq[lq][0]) + (Pyy - nq[lq][1]) * (Pyy - nq[lq][1])) < qRadius) {
                trapped = true;
            }
            
            lq++;
        }
        
        return trapped;
    }        
        
    public void placeCharges(DrawingPanel panel2D) {   
        charge = new InteractiveCircle[q.length];  
        j = 0;
        while(j < q.length) {
            charge[j] = new InteractiveCircle(q[j][0], q[j][1]);
            if(Math.sqrt(q[j][0] * q[j][0] + q[j][1] * q[j][1]) > farthestRadius) {
                farthestRadius = Math.sqrt(q[j][0] * q[j][0] + q[j][1] * q[j][1]);
            }      
                
            panel2D.addDrawable(charge[j]);
            j++;
        }        
        
        panel2D.setPreferredMinMax(-farthestRadius,farthestRadius,-farthestRadius,farthestRadius);
        panel2D.setSquareAspect(true);
        panel2D.repaint();        
    }
    
    /*public void updateChargePosition(double xis, double yp, InteractivePanel panel2D) {
        q[0][0] = xis;
        q[0][1] = yp;
        
        panel2D.repaint();
    }*/
    
    public void createEFieldLines2D_P(javax.swing.JTextField jtext, double precision, int nLines_, DrawingPanel panel2D) {       
        data = new DatasetManager();
        panel2D.addDrawable(data);
        nLines = nLines_;
        Px = new double[nLines];
        Py = new double[nLines];        
        Dx = new double[nLines];
        Dy = new double[nLines];
        Ec = new double[2];
        h = 0;
        while(h < q.length) {
            ang = 0;
            j = 0;
            while(j < nLines) {
                Px[j] = q[h][0] + qRadius * Math.cos(ang);
                Py[j] = q[h][1] + qRadius * Math.sin(ang);
                ang = ang + (2 * Math.PI) / nLines;
                j++;
            }
            
            i = 0;
            while(Math.abs(Px[0]) < farthestRadius*3) {
                data.append(i, Px, Py);
                data.setMarkerShape(i, 6);                
                j = 0;
                while(j < nLines) {
                    Ec = eCalc2D(Px[j], Py[j]);
                    d = Math.sqrt((1 / (precision * precision)) / ((Ec[0] * Ec[0]) + (Ec[1] * Ec[1])));
                    Dx[j] = d * Ec[0];
                    Dy[j] = d * Ec[1];
                    Px[j] = Px[j] + Dx[j];
                    Py[j] = Py[j] + Dy[j];
                    j++;
                }
                
                i++;
            }
            
            h++;
        }
        
        panel2D.repaint();
    }
    
    public void createEFieldLines2D_N(javax.swing.JTextField jtext, double precision, int nLines, DrawingPanel panel2D) {      
        
        //trappedLine = new boolean[nLines];
        data = new DatasetManager();
        panel2D.addDrawable(data);
        Px = new double[nLines];
        Py = new double[nLines];        
        Dx = new double[nLines];
        Dy = new double[nLines];
        Ec = new double[2];
        h = 0;
        //while(h < q.length) {
            //if(q[h][2] > 0) {
            ang = 0;
            j = 0;
            
            while(j < nLines) {
                //trappedLine[j] = false;
                Px[j] = q[h][0] + qRadius * Math.cos(ang);
                Py[j] = q[h][1] + qRadius * Math.sin(ang);
                ang = ang + (2 * Math.PI) / nLines;
                j++;
            }
            
            //allTrapped = false;
            i = 0;
            while(Px[6] < 2.7) {
                data.append(i, Px, Py);
                data.setMarkerShape(i, 6);
                //allTrapped = true;
                j = 0;
                jtext.setText("aqui");
                while(j < nLines) {
                    //if(trappedLine[j] == false) {
                    if(checkEntrapment(Px[j], Py[j]) == false) {
                        //allTrapped = false;
                    Ec = eCalc2D(Px[j], Py[j]);
                    d = Math.sqrt((1 / (precision * precision)) / ((Ec[0] * Ec[0]) + (Ec[1] * Ec[1])));
                    Dx[j] = d * Ec[0];
                    Dy[j] = d * Ec[1];
                    Px[j] = Px[j] + Dx[j];
                    Py[j] = Py[j] + Dy[j];
                    //}
                    }
                    j++;
                }                
                
                i++;
            }
            
            //}
            //h++;
        //}
        
        panel2D.repaint();
        
    }
    
    public void createEFieldVectors2D(boolean real, int grid, double representativeSize, DrawingPanel panel2D) {    
        double Ecc[] = new double[2];
        Data2D data2D_1 = new Data2D(grid, grid, 3);
        data2D_1.setScale(panel2D.getXMin()-0.1*panel2D.getWidth(), panel2D.getXMax()+0.1*panel2D.getWidth(), panel2D.getYMin()-0.1*panel2D.getHeight(), panel2D.getYMax()+0.1*panel2D.getHeight());
        panel2D.setSquareAspect(true);
        double[][][] vector = data2D_1.getData();
       
        for(int iii = 0, nxx = vector.length; iii < nxx; iii++) {
            for(int jjj = 0, nyy = vector[0].length; jjj < nyy; jjj++) {
                double xx = vector[iii][jjj][0];
                double yy = vector[iii][jjj][1];                
                Ecc = eCalc2D(xx,yy);                
                dr = Math.sqrt((representativeSize * representativeSize) / ((Ecc[0] * Ecc[0]) + (Ecc[1] * Ecc[1])));
                if(real == true) {
                    vector[iii][jjj][2] = Math.sqrt(Ecc[0] * Ecc[0] + Ecc[1] * Ecc[1]);
                }
                
                else {
                    vector[iii][jjj][2] = 1;
                }
                    
                vector[iii][jjj][3] = Ecc[0]*dr;
                vector[iii][jjj][4] = Ecc[1]*dr;               
            }
            
        } 
        
        VectorPlot plot1 = new VectorPlot(data2D_1);
        panel2D.addDrawable(plot1);        
        panel2D.repaint();
    }
    
    public void createEFieldIsoPotentialLines2D(DrawingPanel panel2D) {
        Data2D data2D_2 = new Data2D(80,80,1);
        data2D_2.setScale(-11,11,-11,11);//panel2D.getXMin()-0.2*panel2D.getWidth(), panel2D.getXMax()+0.2*panel2D.getWidth(), panel2D.getYMin()-0.2*panel2D.getHeight(), panel2D.getYMax()+0.2*panel2D.getHeight());
        panel2D.setSquareAspect(true);
        double[][][] potential = data2D_2.getData();
        
        for(int ii = 0, nx = potential.length; ii < nx; ii++) {
            for(int jj = 0, ny = potential[0].length; jj < ny; jj++) {
                double x = potential[ii][jj][0];
                double y = potential[ii][jj][1];                
                potential[ii][jj][2] = pCalc2D(x,y);          
            }     
        }        
        
        ContourPlot plot2 = new ContourPlot(data2D_2);
        panel2D.addDrawable(plot2);
        panel2D.repaint();
    }
}
