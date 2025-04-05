import math 
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.style as style
import scipy.special


matplotlib.rcParams['font.family'] = 'DejaVu San'

class Complex_number:  
      
    def __init__(self, real=0, imaginary=0):  
        self.real = real  
        self.imaginary = imaginary  
          
    def __add__(self, other):  
        real = self.real + other.real  
        imag = self.imaginary + other.imaginary  
        return Complex_number(real, imag)   

     # Subtraction
    def __sub__(self, other):
        return Complex_number(self.real - other.real, self.imaginary - other.imaginary) 
  
    def __mul__(self, other):  
        real = self.real * other.real - self.imaginary * other.imaginary  
        imag = self.real * other.imaginary + other.real * self.imaginary  
        return Complex_number(real, imag)    
      
    def __truediv__(self, other):  
        if other.real == 0 and other.imaginary == 0:
            raise ZeroDivisionError("Cannot divide by zero.")
        n = self * Complex_number(other.real, (-1)*other.imaginary)  # Multiply by conjugate
        d = other.real**2 + other.imaginary**2 
        return Complex_number(n.real/d , n.imaginary/d)  
  
    # Exponent
    def __pow__(self, exponent):  
        if exponent == 0:  
            return Complex_number(1, 0) 
        elif exponent < 0:  
            reciprocal = Complex_number(1, 0) / self  # Compute reciprocal
            result = reciprocal ** (-exponent)  # Raise the reciprocal to the positive exponent
        else:  
            r, theta = self.modulus(), self.argument()
            new_r = r**exponent
            new_theta = theta * exponent
            return Complex_number.from_polar(new_r, math.degrees(new_theta)) 

        
    # Conjugate
    def conjugate(self):
        return Complex_number(self.real, -self.imaginary)
    
    # Modulus (Magnitude)
    def modulus(self):
        return math.sqrt(self.real**2 + self.imaginary**2)
    
    # Argument (Angle in radians)
    def argument(self):
        return math.atan2(self.imaginary, self.real)
    
    # Equality check
    def __eq__(self, other):
        return math.isclose(self.real, other.real, rel_tol=1e-9) and math.isclose(self.imaginary, other.imaginary, rel_tol=1e-9)

    
    # Square Root (Returns one root, can be extended for both)
    def sqrt(self):
        r = self.modulus()
        theta = self.argument() / 2
        return Complex_number(math.sqrt(r) * math.cos(theta), math.sqrt(r) * math.sin(theta))
    
    # Rectangular to Polar conversion
    def to_polar(self):
        return (self.modulus(), math.degrees(self.argument()))
    
    # Polar to Rectangular conversion (Static Method)
    @staticmethod
    def from_polar(r, theta):
        theta = math.radians(theta)
        real = r * math.cos(theta)   
        imag = r * math.sin(theta)
        return Complex_number(real, imag)
    
    # Logarithm (Natural Log)
    def log(self):
        r = self.modulus()
        theta = self.argument()
        return Complex_number(math.log(r), theta)
    
    def log_base(self, base):
        return self.log() / math.log(base)

    
    # Complex Trigonometric Functions
    def sin(self):
        return Complex_number(math.sin(self.real) * math.cosh(self.imaginary), math.cos(self.real) * math.sinh(self.imaginary))
    
    def cos(self):
        return Complex_number(math.cos(self.real) * math.cosh(self.imaginary), -math.sin(self.real) * math.sinh(self.imaginary))
    
    def tan(self):
        return self.sin() / self.cos()
    
    # Complex Hyperbolic Functions
    def sinh(self):
        return Complex_number(math.sinh(self.real) * math.cos(self.imaginary), math.cosh(self.real) * math.sin(self.imaginary))
    
    def cosh(self):
        return Complex_number(math.cosh(self.real) * math.cos(self.imaginary), math.sinh(self.real) * math.sin(self.imaginary))
    
    def tanh(self):
        return self.sinh() / self.cosh()
    
    #Inverse trigonometric functions :
    def asin(self):
        i = Complex_number(0, 1)
        one = Complex_number(1, 0)
        z_squared = self * self
        sqrt_expr = (one - z_squared).sqrt()
        inside_log = i * self + sqrt_expr
        result = (inside_log.log()) * Complex_number(0, -1)
        return result

    def acos(self):
        i = Complex_number(0, 1)
        one = Complex_number(1, 0)
        z_squared = self * self
        sqrt_expr = (z_squared - one).sqrt()
        inside_log = self + sqrt_expr
        result = (inside_log.log()) * Complex_number(0, -1)
        return result

    def atan(self):
        i = Complex_number(0, 1)
        one = Complex_number(1, 0)
        numerator = (one - i * self).log()
        denominator = (one + i * self).log()
        result = (numerator - denominator) * Complex_number(0, 0.5)
        return result

    
    # Complex n-th Roots (Returns all n roots)
    def nth_roots(self, n):
        roots = []
        r = self.modulus() ** (1/n)
        theta = self.argument()
        for k in range(n):
            angle = (theta + 2 * math.pi * k) / n
            roots.append(Complex_number(r * math.cos(angle), r * math.sin(angle)))
        return roots

    # Exponential function (e^z) : 
    def exp(self):
        real_part = math.exp(self.real) * math.cos(self.imaginary)
        imag_part = math.exp(self.real) * math.sin(self.imaginary)
        return Complex_number(real_part, imag_part)
    
    # Gamma Function for Complex Numbers (Extension of factorial)
    def gamma(self):
        real_part = scipy.special.gamma(self.real) * math.cos(self.imaginary)
        imag_part = scipy.special.gamma(self.real) * math.sin(self.imaginary)
        return Complex_number(real_part, imag_part)


    def __round__(self, ndigits=10):
        return Complex_number(round(self.real, ndigits), round(self.imaginary, ndigits))

    def __ceil__(self):
        return Complex_number(math.ceil(self.real), math.ceil(self.imaginary))

    def __floor__(self):    
        return Complex_number(math.floor(self.real), math.floor(self.imaginary))


  
    def __repr__(self):  
        # For better printing of complex numbers  
        if self.real == 0 and self.imaginary == 0:
            return "0"
        elif self.real == 0:
            return f"{self.imaginary}i"
        elif self.imaginary == 0:
            return f"{self.real}"
        elif self.imaginary > 0:
            return f"{self.real} + {self.imaginary}i"
        else:
            return f"{self.real} - {-self.imaginary}i" 


    # Required VISUALIZATION Functions:

    def plot(self, others=None, title="Complex Numbers in the Complex Plane", show=True):
        """
        Plot this complex number or a list of complex numbers in the complex plane with a professional design.
        If 'others' is provided, it should be a list of Complex_number instances to plot alongside self.
        """
        # Create a figure with a specific size and DPI for clarity
        plt.figure(figsize=(8, 8), dpi=100)

        # Set a light background and professional style
        plt.style.use('seaborn-v0_8-whitegrid')

        # Plot axes with a subtle gray color
        plt.axhline(0, color='gray', lw=1, alpha=0.5)
        plt.axvline(0, color='gray', lw=1, alpha=0.5)

        # Colors for multiple complex numbers (cycler-like behavior)
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']  # Blue, Orange, Green, Red, Purple

        # If 'others' is None, plot only self; otherwise, include all in the list
        if others is None:
            numbers = [self]
        else:
            if not isinstance(others, list):
                raise ValueError("'others' must be a list of Complex_number instances")
            if not all(isinstance(z, Complex_number) for z in others):
                raise ValueError("All elements in 'others' must be Complex_number instances")
            numbers = [self] + others

        # Plot each complex number
        for i, z in enumerate(numbers):
            color = colors[i % len(colors)]  # Cycle through colors if more than 5 numbers
            # Plot point
            plt.plot(z.real, z.imaginary, 'o', color=color, markersize=10, 
                    markeredgecolor='black', markeredgewidth=1, label=f'{z}')
            # Plot vector
            plt.quiver(0, 0, z.real, z.imaginary, angles='xy', scale_units='xy', scale=1, 
                    color=color, alpha=0.8, width=0.005, headwidth=5)

            # Get polar form
            r, theta = z.to_polar()

            # Annotate with both rectangular and polar form
            annotation_text = f'{z.real:.2f} {"+" if z.imaginary >= 0 else "-"} {abs(z.imaginary):.2f}i\n' \
                            f'r = {r:.2f}, θ = {theta:.1f}°'

            plt.annotate(annotation_text,
                        xy=(z.real, z.imaginary),
                        xytext=(10, 10), textcoords='offset points',
                        fontsize=10, color=color,
                        bbox=dict(boxstyle="round,pad=0.3", edgecolor="none", facecolor="white", alpha=0.8))

        # Dynamic axis limits based on the maximum modulus
        max_modulus = max(z.modulus() for z in numbers)
        max_val = max(1.5 * max_modulus, 1)  # Ensure a minimum range for small numbers
        plt.xlim(-max_val, max_val)
        plt.ylim(-max_val, max_val)

        # Dynamic tick marks (adjust based on maximum magnitude)
        tick_step = max_val / 4  # Roughly 4 ticks on each side of zero
        if tick_step < 0.1:
            tick_step = round(tick_step, 2)  # Small numbers get finer ticks
        else:
            tick_step = round(tick_step, 1)  # Larger numbers get coarser ticks
        x_ticks = [i * tick_step for i in range(int(-max_val / tick_step), int(max_val / tick_step) + 1)]
        y_ticks = x_ticks  # Symmetric for real and imaginary axes
        plt.xticks(x_ticks, fontsize=10)
        plt.yticks(y_ticks, fontsize=10)

        # Customize labels and title with a professional font
        plt.xlabel('Real Axis', fontsize=12, weight='bold')
        plt.ylabel('Imaginary Axis', fontsize=12, weight='bold')
        plt.title(title, fontsize=14, weight='bold', pad=15)

        # Add a legend with a clean design
        plt.legend(loc='upper left', fontsize=10, frameon=True, edgecolor='black', facecolor='white', framealpha=0.9)

        # Ensure equal aspect ratio for accurate representation
        plt.axis('equal')

        # Tight layout for a polished look
        plt.tight_layout()

        if show:
            plt.show()
        else:
            return plt  # Allow further customization if not shown

    
    def plot_polar(self, others=None, title="Complex Numbers in Polar Form", show=True):
        """
        Visualize this complex number or a list of complex numbers in polar form, showing magnitude and angle.
        If 'others' is provided, it should be a list of Complex_number instances to plot alongside self.
        """
        plt.figure(figsize=(8, 8), dpi=100)
        plt.style.use('seaborn-v0_8-whitegrid')
        
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']  # Blue, Orange, Green, Red, Purple
        
        if others is None:
            numbers = [self]
        else:
            if not isinstance(others, list) or not all(isinstance(z, Complex_number) for z in others):
                raise ValueError("'others' must be a list of Complex_number instances")
            numbers = [self] + others
        
        for i, z in enumerate(numbers):
            color = colors[i % len(colors)]
            r, theta = z.to_polar()
            theta_rad = math.radians(theta)
            
            # Plot the line from origin to point
            plt.polar([0, theta_rad], [0, r], '-', color=color, lw=2, alpha=0.8)
            # Plot the point
            plt.plot(theta_rad, r, 'o', color=color, markersize=10, markeredgecolor='black', 
                    markeredgewidth=1, label=f'{z} (r={r:.2f}, θ={theta:.2f}°)')
            
            # Annotate near the point (convert to Cartesian for positioning)
            x, y = r * math.cos(theta_rad), r * math.sin(theta_rad)
            plt.annotate(f'r = {r:.2f}\nθ = {theta:.1f}°', 
                        xy=(theta_rad, r), xycoords='data', 
                        xytext=(10, 10), textcoords='offset points', 
                        fontsize=10, color=color, 
                        bbox=dict(boxstyle="round,pad=0.3", edgecolor="none", facecolor="white", alpha=0.8))
        
        plt.title(title, fontsize=14, weight='bold', pad=15)
        plt.legend(loc='upper left', fontsize=10, frameon=True, edgecolor='black', facecolor='white', framealpha=0.9)
        plt.tight_layout()
        
        if show:
            plt.show()
        else:
            return plt
        

    def plot_nth_roots(self, n, title="Nth Roots of Complex Number", show=True):
        """
        Plot all nth roots of this complex number in the complex plane with a professional design.
        """
        plt.figure(figsize=(8, 8), dpi=100)
        plt.style.use('seaborn-v0_8-whitegrid')
        
        plt.axhline(0, color='gray', lw=1, alpha=0.5)
        plt.axvline(0, color='gray', lw=1, alpha=0.5)
        
        roots = self.nth_roots(n)
        color = '#1f77b4'  # Single color for roots (blue)
        
        for i, root in enumerate(roots):
            # Plot point
            plt.plot(root.real, root.imaginary, 'o', color=color, markersize=8, 
                    markeredgecolor='black', markeredgewidth=1, label=f'Root {i}: {root}' if i == 0 else None)
            # Plot vector
            plt.quiver(0, 0, root.real, root.imaginary, angles='xy', scale_units='xy', scale=1, 
                    color=color, alpha=0.6, width=0.004, headwidth=5)
            # Annotate
            r, theta = root.to_polar()
            plt.annotate(f'r = {r:.2f}\nθ = {theta:.1f}°', 
                        xy=(root.real, root.imaginary), 
                        xytext=(10, 10), textcoords='offset points', 
                        fontsize=9, color=color, 
                        bbox=dict(boxstyle="round,pad=0.3", edgecolor="none", facecolor="white", alpha=0.8))
        
        # Plot the original number for reference
        plt.plot(self.real, self.imaginary, 'o', color='#ff7f0e', markersize=10, 
                markeredgecolor='black', markeredgewidth=1, label=f'Original: {self}')
        plt.quiver(0, 0, self.real, self.imaginary, angles='xy', scale_units='xy', scale=1, 
                color='#ff7f0e', alpha=0.8, width=0.005, headwidth=5)
        
        # Dynamic axis limits
        max_modulus = max(max(root.modulus() for root in roots), self.modulus())
        max_val = max(1.5 * max_modulus, 1)
        plt.xlim(-max_val, max_val)
        plt.ylim(-max_val, max_val)
        
        # Dynamic ticks
        tick_step = max_val / 4
        if tick_step < 0.1:
            tick_step = round(tick_step, 2)
        else:
            tick_step = round(tick_step, 1)
        x_ticks = [i * tick_step for i in range(int(-max_val / tick_step), int(max_val / tick_step) + 1)]
        plt.xticks(x_ticks, fontsize=10)
        plt.yticks(x_ticks, fontsize=10)
        
        plt.xlabel('Real Axis', fontsize=12, weight='bold')
        plt.ylabel('Imaginary Axis', fontsize=12, weight='bold')
        plt.title(title, fontsize=14, weight='bold', pad=15)
        plt.legend(loc='upper left', fontsize=10, frameon=True, edgecolor='black', facecolor='white', framealpha=0.9)
        plt.axis('equal')
        plt.tight_layout()
        
        if show:
            plt.show()
        else:
            return plt
        

    def plot_with(self, other, operation='add', title="Operation Result", show=True):
        """
        Plot this complex number, another complex number, and their operation result.
        Supported operations: 'add', 'sub', 'mul', 'div'.
        'other' can be a single Complex_number or a list for multiple operations.
        """
        plt.figure(figsize=(8, 8), dpi=100)
        plt.style.use('seaborn-v0_8-whitegrid')
        
        plt.axhline(0, color='gray', lw=1, alpha=0.5)
        plt.axvline(0, color='gray', lw=1, alpha=0.5)
        
        colors = {'self': '#1f77b4', 'other': '#ff7f0e', 'result': '#2ca02c'}  # Blue, Orange, Green
        
        if isinstance(other, Complex_number):
            others = [other]
        elif isinstance(other, list) and all(isinstance(z, Complex_number) for z in other):
            others = other
        else:
            raise ValueError("'other' must be a Complex_number or a list of Complex_number instances")
        
        # Plot self
        plt.plot(self.real, self.imaginary, 'o', color=colors['self'], markersize=10, 
                markeredgecolor='black', markeredgewidth=1, label=f'Self: {self}')
        plt.quiver(0, 0, self.real, self.imaginary, angles='xy', scale_units='xy', scale=1, 
                color=colors['self'], alpha=0.8, width=0.005, headwidth=5)
        
        
        z = self
        color = colors['self']
        annotation = (
            f'{z.real:.2f} {"+" if z.imaginary >= 0 else "-"} {abs(z.imaginary):.2f}i\n'
            f'r = {z.modulus():.2f}, θ = {z.to_polar()[1]:.1f}°'
        )
        plt.annotate(
            annotation,
            xy=(z.real, z.imaginary),
            xytext=(10, 10),
            textcoords='offset points',
            fontsize=10,
            color=color,
            bbox=dict(boxstyle="round,pad=0.3", edgecolor="none", facecolor="white", alpha=0.85)
        )
        
        results = []
        for i, o in enumerate(others):
            color_other = colors['other'] if i == 0 else '#d62728'
            plt.plot(o.real, o.imaginary, 'o', color=color_other, markersize=10, 
                    markeredgecolor='black', markeredgewidth=1, label=f'Other {i+1}: {o}')
            plt.quiver(0, 0, o.real, o.imaginary, angles='xy', scale_units='xy', scale=1, 
                    color=color_other, alpha=0.8, width=0.005, headwidth=5)
            
           
            z = o
            color = color_other
            annotation = (
                f'{z.real:.2f} {"+" if z.imaginary >= 0 else "-"} {abs(z.imaginary):.2f}i\n'
                f'r = {z.modulus():.2f}, θ = {z.to_polar()[1]:.1f}°'
            )
            plt.annotate(
                annotation,
                xy=(z.real, z.imaginary),
                xytext=(10, 10),
                textcoords='offset points',
                fontsize=10,
                color=color,
                bbox=dict(boxstyle="round,pad=0.3", edgecolor="none", facecolor="white", alpha=0.85)
            )
            
            # Compute result
            if operation == 'add':
                result = self + o
            elif operation == 'sub':
                result = self - o
            elif operation == 'mul':
                result = self * o
            elif operation == 'div':
                result = self / o
            else:
                raise ValueError("Unsupported operation. Use 'add', 'sub', 'mul', or 'div'")
            
            results.append(result)
            plt.plot(result.real, result.imaginary, 'o', color=colors['result'], markersize=10, 
                    markeredgecolor='black', markeredgewidth=1, label=f'Result {i+1}: {result}' if i == 0 else None)
            plt.quiver(0, 0, result.real, result.imaginary, angles='xy', scale_units='xy', scale=1, 
                    color=colors['result'], alpha=0.8, width=0.005, headwidth=5)
            
           
            z = result
            color = colors['result']
            annotation = (
                f'{z.real:.2f} {"+" if z.imaginary >= 0 else "-"} {abs(z.imaginary):.2f}i\n'
                f'r = {z.modulus():.2f}, θ = {z.to_polar()[1]:.1f}°'
            )
            plt.annotate(
                annotation,
                xy=(z.real, z.imaginary),
                xytext=(10, 10),
                textcoords='offset points',
                fontsize=10,
                color=color,
                bbox=dict(boxstyle="round,pad=0.3", edgecolor="none", facecolor="white", alpha=0.85)
            )
        
        max_modulus = max(max(self.modulus(), o.modulus(), r.modulus()) for o, r in zip(others, results))
        max_val = max(1.5 * max_modulus, 1)
        plt.xlim(-max_val, max_val)
        plt.ylim(-max_val, max_val)
        
        tick_step = max_val / 4
        tick_step = round(tick_step, 2) if tick_step < 0.1 else round(tick_step, 1)
        x_ticks = [i * tick_step for i in range(int(-max_val / tick_step), int(max_val / tick_step) + 1)]
        plt.xticks(x_ticks, fontsize=10)
        plt.yticks(x_ticks, fontsize=10)
        
        plt.xlabel('Real Axis', fontsize=12, weight='bold')
        plt.ylabel('Imaginary Axis', fontsize=12, weight='bold')
        plt.title(title, fontsize=14, weight='bold', pad=15)
        plt.legend(loc='upper left', fontsize=10, frameon=True, edgecolor='black', facecolor='white', framealpha=0.9)
        plt.axis('equal')
        plt.tight_layout()
        
        if show:
            plt.show()
        else:
            return plt


# TESTING:

if __name__ == "__main__":
    z1 = Complex_number(2, 7)
    z2 = Complex_number(1, -3)
    z3 = Complex_number(-4, 2)
    
    # Test plot_polar
    z2.plot();
    z1.plot_polar(title="Single Polar: z1")
    z1.plot_polar([z2, z3], title="Multiple Polar: z1, z2, z3")
    
    # Test plot_nth_roots
    z1.plot_nth_roots(3, title="Cube Roots of z1")
    
    # Test plot_with
    z1.plot_with(z2, operation='add', title="z1 + z2")
    z1.plot_with([z2, z3], operation='mul', title="z1 * [z2, z3]")