# AAE6102_Assignment

This respository provides the calculation of obtaining the GPS user position from receiver and ephermeris data. The aim of this assignment is to research and apply the knowledge learnt in AAE6102.

# To Use
1. Please clone or download this respository
2. Open matlab and add the parent folder to the matlab path
3. Run the main.m in matlab

# Methodology
The following flowchart shows the data relation and the algorithm for estimating the position of the GPS receiver. The first step is using receiver and ephemeris data to estimate the position of the satellite and the clock offset. This can be subdivided into three sections, the approximate satellite position is calculated using orbit parameters, the earth rotation is applied to correct the satellite positions, the broadcast satellite clock error is calculated based on the eccentric anomaly. The second step is to determine the receiver position based on the satellite position and clock offset. Firstly, the tropospheric delay is calculated based on standard atmospheric conditions using the Hopfield model. Then the receiver position calculation is linearized with the pseudo-range observation and calculation, the receiver 3D position and clock offset estimation is solved through the least-squares iteration. Finally, the receiver position and clock offset is calculated.
![image](https://user-images.githubusercontent.com/65110263/139574871-066e2f0e-3dda-45db-a456-c71f34571dca.png)

# Results
Initial position in earth-centred-earth-fixed (ECEF) coordinates (in meter) is **[-2694685.473, -4293642.366, 3857878.924]**.

With the satellite coordinates, measured pseudo-ranges, and the estimated receiver position, the most probable position of the receiver can now be found. In a least squares adjustment (Hofmann-Wellenhof et al. 2001, p. 256), the computed ranges, observed ranges, and the receiver clock bias are applied to a set of linearized model equations (equal to the number of satellites involved) to find the final receiver position. This solution is an approximation, so an iterative approach is required whereupon the most recent solution becomes the initial guess, and the process is repeated until the desired accuracy is obtained. If more than four satellites are being observed, the problem is over-determined and can be solved in a least squares sense to yield an optimal receiver location.

After several iteration, the ECEF solution (in meter) is **[-2700419.67455933	-4292538.83987527	3855266.39572346]**. In WGS84 coordinates format with latitude (degree), longitude (degree), and altitude(meter), it is **[37.4281°, -122.1738°, 60.5360m]**. The estimated receiver clock offset is **519468.902m (0.00173276s)**. The estimated position error is **29.11752m**.

If the code excuted successfully, following figures will be shown.


![image](https://user-images.githubusercontent.com/65110263/139575127-83812d5d-a829-41d1-bd4d-612381ce48b2.png) ![image](https://user-images.githubusercontent.com/65110263/139575131-37c9b956-782c-450d-9771-9fd40402222a.png) ![image](https://user-images.githubusercontent.com/65110263/139575133-ba711a33-3869-4783-b7a1-afeb1802f5ef.png)

# Author
For any issues, please contact Max via email 21106098R@connect.polyu.hk
