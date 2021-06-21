/*
 * EField3D.java
 *
 * Created on 13 de Junho de 2004, 15:40
 */

/**
 *
 * @author  user
 */

//import org.opensourcephysics.displayejs.*;
//import org.opensourcephysics.ejs.display.*;
import org.opensourcephysics.displayejs.*;
public class EField3D {
    private int g, h, i, j, m, YLines = 6, ZLines = 6;
    private double q[][], x[][], y[][], z[][];
    private double Ec[], E[], Ex[], Ey[], Ez[], r[], ExEyEzSum[];
    private double dx, dy, dz, d, precision = 0.05, angY, angZ, size, ExSum, EySum, EzSum, pX, pY, pZ, k = 8.987551788 * Math.pow(10, 9);
    private DrawingPanel3D panel3D = new DrawingPanel3D();
    private InteractiveParticle charge[], dot[][][];
    
    public EField3D(double[][] q_) {
        q = new double[q_.length+1][4];
        q[0][0] = 0;
        q[0][1] = 0;
        q[0][2] = 0;
        q[0][3] = 1;
        g = 0;
        while(g < q_.length) {
            q[g+1][0] = q_[g][0];
            q[g+1][1] = q_[g][1];
            q[g+1][2] = q_[g][2];
            q[g+1][3] = q_[g][3];
            g++;
        }        
               
        eCalc3D(pX, pY, pZ);        
        placeCharges(panel3D); 
        createEFieldLines3D_P(panel3D);        
    }
    
    public double[] eCalc3D(double pX, double pY, double pZ) {
        r = new double[q.length];
        E = new double[q.length];
        Ex = new double[q.length];
        Ey = new double[q.length];
        Ez = new double[q.length];
        ExEyEzSum = new double[3];
        
        m = 0;
        ExSum = 0;
        EySum = 0;
        EzSum = 0;
        while(m < q.length) {
            r[m] = Math.sqrt((pX - q[m][0]) * (pX - q[m][0]) + (pY - q[m][1]) * (pY - q[m][1]) + (pZ - q[m][1]) * (pZ - q[m][1]));
            E[m] = (k * q[m][2]) / (r[m] * r[m]);
            Ex[m] = E[m] * ((pX - q[m][0]) / r[m]);
            Ey[m] = E[m] * ((pY - q[m][1]) / r[m]);
            Ez[m] = E[m] * ((pZ - q[m][1]) / r[m]);
            ExSum = ExSum + Ex[m];
            EySum = EySum + Ey[m];
            EzSum = EzSum + Ez[m];
            m++;
        }
        
        ExEyEzSum[0] = ExSum;
        ExEyEzSum[1] = EySum;
        ExEyEzSum[2] = EzSum;
        return ExEyEzSum;
    }
    
    public void placeCharges(DrawingPanel3D panel3D) {
        charge = new InteractiveParticle[q.length];
        
        g = 0;
        while(g < q.length) {
            charge[g] = new InteractiveParticle();
            g++;
        }
        
        j = 0;
        while(j < q.length) {
            charge[j].setXYZ(q[j][0], q[j][1], q[j][2]);
            charge[j].setSizeX(0.3);
            charge[j].setSizeY(0.3);
            charge[j].setSizeZ(0.3);           
            panel3D.addDrawable(charge[j]);            
            j++;
        }  
        panel3D.setSquareAspect(true);
        panel3D.repaint();
        
    }
    
    public void createEFieldLines3D_P(/*javax.swing.JTextField jtext, double precision, int nLines, */DrawingPanel3D panel3D) {      
        x = new double[YLines][ZLines];
        y = new double[YLines][ZLines];
        z = new double[YLines][ZLines];
        dot = new InteractiveParticle[120][YLines][ZLines];
        Ec = new double[3];
        h = 0;
        //while(h < q.length) {
            size = charge[h].getSizeX();
            angZ = 0;
            i = 0;
            j = 0;
            while(j < ZLines) {
                angY = 0;
                g = 0;
                while(g < YLines) {                    
                    x[j][g] = size * Math.cos(angZ) * Math.cos(angY);
                    y[j][g] = size * Math.cos(angZ) * Math.sin(angY);
                    z[j][g] = size * Math.sin(angZ);
                    dot[i][j][g] = new InteractiveParticle();                    
                    dot[i][j][g].setXYZ(x[j][g], y[j][g], z[j][g]);
                    dot[i][j][g].setShapeType(0);
                    panel3D.addDrawable(dot[i][j][g]);
                    angY = angY + (2 * Math.PI) / YLines;
                    g++;
                }
                
                angZ = angZ + (2 * Math.PI) / ZLines;
                j++;
            }            
            i = 0;
            while(x[0][0] < 5) {
            angZ = 0;
            j = 0;
            while(j < ZLines) {
                angY = 0;
                g = 0;
                while(g < YLines) {                    
                    Ec = eCalc3D(x[j][g], y[j][g], z[j][g]);
                    d = Math.sqrt((1 / (precision * precision)) / (Ec[0] * Ec[0] + Ec[1] * Ec[1] + Ec[2] * Ec[2]));                    
                    dx = d * Ec[0];
                    dy = d * Ec[1];
                    dz = d * Ec[2];
                    x[j][g] = x[j][g] + dx;
                    y[j][g] = y[j][g] + dy;
                    z[j][g] = z[j][g] + dz;
                    dot[i][j][g] = new InteractiveParticle();
                    dot[i][j][g].setShapeType(0);
                    dot[i][j][g].setXYZ(x[j][g], y[j][g], z[j][g]);
                    panel3D.addDrawable(dot[i][j][g]);
                    angY = angY + (2 * Math.PI) / YLines;
                    g++;
                }
                
                angZ = angZ + (2 * Math.PI) / ZLines;
                j++;
            }
            i++;
            }
            //panel3D.repaint();
    }    
}
