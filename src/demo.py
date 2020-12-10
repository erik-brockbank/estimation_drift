import pygame
import math, sys, os, shutil
from experimenter import *
from numpy import random as rand
import random
import time

white = (255,255,255)
black = (0,0,0)

diameter = 20
n = 100 # tweak this for modified demo
#calibrations = 3

def waitForKey():
	pressed = False
	while(not pressed):
		for event in pygame.event.get():
			if event.type == KEYDOWN:
				if event.key == K_SPACE:
					pressed = 1
				elif event.key == K_ESCAPE:
					quit('I pressed ESCAPE!')
	return(pressed)

def buildConfiguration(n):
    configuration = []
    for i in range(0,n):
        x = random.randint(int(.05*dispSize[0]), dispSize[0] - int(.05*dispSize[0]) - diameter)
        y = random.randint(int(.05*dispSize[1]), dispSize[1] - int(.05*dispSize[1]) - diameter)
        rect = pygame.Rect(x,y,diameter,diameter)
        if rect.collidelist(configuration) == -1:
            configuration.append(rect)
    return configuration

def drawDots(drawer, configuration):
    # makes dots
    drawer.screenFill((white))
    for rect in configuration:
        # pygame.draw.circle(drawer.screen, (200, 0, 0), rect.center, diameter/2, 0)
        pygame.draw.circle(drawer.screen, (200, 0, 0), rect.center, int(diameter/2), 0) #eb
    drawer.screenFlip()
    waitForKey()

# Initialize the game engine
pygame.init()
# Get subject ID and set up a file
# Initialize drawer
drawer = Experimenter()

dispSize = [int(drawer.screen.get_width()*.95), int(drawer.screen.get_height()*.95)]
# start timing the experiment

drawDots(drawer, buildConfiguration(n))

pygame.quit()
sys.exit(0)
