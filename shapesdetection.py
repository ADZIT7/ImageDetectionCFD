import cv2
import numpy as np

# Load image
image = cv2.imread('circle.png')
gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)


# Edge detection
edges = cv2.Canny(gray, threshold1=50, threshold2=150)


# Find contours
contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)


# Draw and extract (x, y) points
for contour in contours:
    for point in contour:
        x, y = point[0]
        print(f'{x} {y}')
        cv2.circle(image, (x, y), 2, (0, 255, 0), -1)  # optional: draw points


# Show result
cv2.imshow("Detected Points", image)
cv2.waitKey(0)
cv2.destroyAllWindows()
