import csv
import matplotlib.pyplot as plt

def plot_histogram(csv_file, output_image):
    avg_distances = []

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            avg_distances.append(float(row['Avg_Distance']))

    plt.figure(figsize=(10, 6))
    plt.hist(avg_distances, bins=20, color='blue', edgecolor='black')
    plt.title('Histogram of Average Distances')
    plt.xlabel('Average Distance (Å)')
    plt.ylabel('Frequency')
    plt.grid(True)
    
    # Save the plot
    plt.savefig(output_image)

def plot_boxplot(csv_file, output_image):
    avg_distances = []

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            avg_distances.append(float(row['Avg_Distance']))

    plt.figure(figsize=(10, 6))
    plt.boxplot(avg_distances, vert=False)
    plt.title('Box Plot of Average Distances')
    plt.xlabel('Average Distance (Å)')
    plt.grid(True)
    
    # Save the plot
    plt.savefig(output_image)

def plot_bar(csv_file, output_image):
    num_aa = []

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            num_aa.append(int(row['Num_AA']))

    plt.figure(figsize=(10, 6))
    plt.hist(num_aa, bins=range(min(num_aa), max(num_aa) + 2), color='blue', edgecolor='black', align='left')
    plt.title('Bar Plot of Number of Amino Acids in Groups')
    plt.xlabel('Number of Amino Acids')
    plt.ylabel('Frequency')
    plt.grid(True)
    
    # Save the plot
    plt.savefig(output_image)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Plot the output of the graph script')
    parser.add_argument('csv_file', help='CSV file containing the output data')
    parser.add_argument('output_image_histogram', help='Output image file to save the histogram plot')
    parser.add_argument('output_image_boxplot', help='Output image file to save the box plot')
    parser.add_argument('output_image_bar', help='Output image file to save the bar plot')
    args = parser.parse_args()

    plot_histogram(args.csv_file, args.output_image_histogram)
    plot_boxplot(args.csv_file, args.output_image_boxplot)
    plot_bar(args.csv_file, args.output_image_bar)