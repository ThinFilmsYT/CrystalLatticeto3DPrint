import numpy as np
import random
import matplotlib.pyplot as plt

def spherical_to_cartesian_vector(azimuth, polar):
    """Converts spherical angles to a cartesian unit vector."""
    x = np.sin(polar) * np.cos(azimuth)
    y = np.sin(polar) * np.sin(azimuth)
    z = np.cos(polar)
    return np.array([x, y, z])

def setup_lines(angle_list):
    """Takes a list of (azimuth, polar) tuples and returns a list of unit vectors."""
    return [spherical_to_cartesian_vector(az, pol) for pol, az in angle_list]

# --- Objective Function ---
def min_angle_to_plane(normal, lines):
    """Computes the minimal angle between the plane (normal) and all lines.
        Assumes lines are one-directional, so does not take absolute value.
    """
    min_angle = np.pi #initialize min angle
    normal /= np.linalg.norm(normal)
    for line in lines:
        line /= np.linalg.norm(line)
        cos_theta = abs(np.dot(normal, line)) # normal dot vector = cos(angle between them) (wlog define positive)
        angle_to_plane = np.pi / 2 - np.arccos(cos_theta) # normal is 90 degrees from plane, so this calculates angle to plane
        if angle_to_plane < min_angle:
            min_angle = angle_to_plane
    return min_angle

# --- Simulated Annealing ---
def simulated_annealing(lines, steps=10000, initial_temp=1.0, cooling_rate=0.999,intial_perturb = 1):
    #initialize useful functions
    def random_direction():
        vec = np.random.normal(size=3)
        return vec / np.linalg.norm(vec)

    def perturb(normal,scale = 0.1):
        perturbation = np.random.normal(scale=scale, size=3)
        new_normal = normal + perturbation
        return new_normal / np.linalg.norm(new_normal)
    
    # initialize random guess for starting plane
    current = random_direction()
    current_score = min_angle_to_plane(current, lines)
    best = current
    best_score = current_score
    T = initial_temp
    perturbation = intial_perturb

    for _ in range(steps):
        candidate = perturb(current,perturbation)
        candidate_score = min_angle_to_plane(candidate, lines)
        # Always accept candidate solution if candidate is better.
        # Alternately, accept candidate solution with some randomness.
        # Just as in a real physical system, probabilty of a state is proportional
        # to the Boltzmann factor exp(energy of state / T) #k_B = 1
        if candidate_score > current_score or random.random() < np.exp((candidate_score - current_score) / T): 
            current = candidate
            current_score = candidate_score
            # if score changed check if it is best so far
            if current_score > best_score:
                best = current
                best_score = current_score

        T *= cooling_rate #exponentially decrease temp as in Newton's Law of Cooling
        perturbation *= cooling_rate

    return best, best_score

def render_plane_and_lines(spherical_lines, plane_normal):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    lines = setup_lines(spherical_lines)
    # Plot each line (as a ray from origin)
    for line in lines:
        ax.plot([0, line[0]], [0, line[1]], [0, line[2]], 'r')

    # Plot the plane (through origin, oriented by normal)
    plane_normal = plane_normal / np.linalg.norm(plane_normal)
    d = 1.5  # extent of the plane
    xx, yy = np.meshgrid(np.linspace(-d, d, 10), np.linspace(-d, d, 10))

    # Solve for zz in the plane equation n . (x,y,z) = 0
    a, b, c = plane_normal
    zz = (-a * xx - b * yy) / c if c != 0 else 0 * xx  # avoid divide-by-zero

    ax.plot_surface(xx, yy, zz, alpha=0.3, color='blue')

    # Style plot
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_zlim([-1.5, 1.5])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Optimal Plane and Directional Lines')
    plt.tight_layout()
    plt.show()

def simulated_anneal_atoms(bonds):
    angle_list = []
    for bond in bonds:
        angle_list.append((bond[1],bond[2]))
    lines = setup_lines(angle_list)
    return simulated_annealing(lines)
    



if __name__ == "__main__":
    # Define lines by (polar, azimuth) angles in radians
    angle_list = [
        (np.pi / 4,0),
        (np.pi / 4,np.pi / 2),
        (np.pi / 4,np.pi),
        (np.pi / 4,3*np.pi / 2),
        # (np.pi / 2,0),
    ]

    # angle_list = [#tetragonally coordinated
    #     (.5*np.arccos(-1/3),0),
    #     (.5*np.arccos(-1/3),np.pi),
    #     (np.pi - .5 * np.arccos(-1/3),np.pi / 2),
    #     (np.pi - .5 * np.arccos(-1/3),3*np.pi / 2),
    # ]

    angle_list = [#octahedrally coordinated
        (0, 0), #polar,azimuth
        (np.pi, 0),
        (np.pi/2, 0),
        (np.pi/2, np.pi/2),
        (np.pi/2, np.pi),
        (np.pi/2, 3*np.pi/2),
    ]

    lines = setup_lines(angle_list)
    optimal_plane_normal, optimal_angle = simulated_annealing(lines)
    
    print("Optimal plane normal vector:", optimal_plane_normal)
    print("Maximum minimal angle to any line (in degrees):", np.degrees(optimal_angle))
    print(optimal_angle/(.5*(np.pi-np.arccos(-1/3))))
    render_plane_and_lines(angle_list,optimal_plane_normal)
