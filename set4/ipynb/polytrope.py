def rk4(n,h,disp = False,ndisp=100):
    """
    Follow the notation and recipe from HWT chapter 7:
        y = theta_n
        x = xi
        z = dy/dx
        
    """
    
    _x = []
    _y = []
    _z = []
    
    
    # Start the integration just off the origin and get y and y' from the taylor
    # expansion for y
    
    x0 = .01367 * h   
    y = 1. - (x0**2)/6. + (n * x0**4 / 120.) - (n*(8*n -5)*x0**6/15120.)
    z = -(x0/3.) + (n*x0**3/30.) - ( 6*n*(8*n -5)*x0**5/15120. )
    
    if disp: print 'z0 = %f, y0 = %f' % (z,y)
       
    i = 0
    while y > 0:
        
        try:
        
            x1 =  x0 + h*i
            x23 = x1 + .5*h
            x4 = x1 + h
        
            k1 = h*z
            l1 = -h*( (y**n) + (2./x1)*z)
        
            z2 = z + (.5 * l1)
            k2 = h * z2
            l2 = -h * ( (y + .5 * k1)**n  + (2./x23) * z2  )
        
            z3 = z + (.5*l2)
            k3 = h * z3
            l3 = -h * (  (y + .5 * k2)**n  + (2./x23) * z3 )
        
            z4 = z + (.5*l3)
            k4 = h * (z + l3)
            l4 = -h * ( (y + k3)**n  + (2./x4) * z4 )
        
            dz = (l1 + 2*l2 + 2*l3 + l4) / 6.
            dy = (k1 + 2*k2 + 2*k3 + k4) / 6.        
        
        
            if (not(i%ndisp)) and disp:
                print'--------'
                print 'x1',x1
                print 'l',l1,l2,l3,l4
                print 'k',k1,k2,k3,k4
                print 'dy,dz',dy,dz
                print 'y,z',y,z
        
            z += dz
            y += dy
        
            i += 1
        
            _y.append(y)
            _x.append(x1)
            _z.append(z)    
    
            if disp: print 'x_1 = %f, y_1 = %f, z_1 = %f, D = %f' % (_x[-1],_y[-1],_z[-1], -_x[-1]**2*_z[-1])
            
        except:
            break
    
    return np.array(_x), np.array(_y), np.array(_z)
