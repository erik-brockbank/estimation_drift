from __future__ import division
# import OpenGL.GL as gl # eb commented out, doesn't seem to be used and breaks otherwise
import pygame
from pygame.locals import *
import pygame.image
from numpy import *
import colorsys

pygame.font.init()


class Experimenter(object):
    def __init__(self):
        self.screen = pygame.display.set_mode((1024,768))#16) # pygame.FULLSCREEN
        # pygame.mouse.set_visible(False)
        self.center = (512,384)
        self.deg2pix = 30
        self.xlim = (0, 1024)
        self.ylim = (0, 768)
        self.font = pygame.font.Font(None, 40)
        self.textDic = {}

    def writeInstructions(self, text, offset=-450):
        counter = 0
        for line in text:
            # if not self.textDic.has_key(line):
            if not (line in self.textDic): #eb
                self.textDic[line] = self.font.render(line, 1, (0,0,0))
            self.screen.blit(self.textDic[line], (self.center[0]+offset, self.center[1]-350+counter*50))
            counter += 1

    def writeScore(self, newpoints, newtotal):
        self.screenFill()
        txtFPS = self.font.render('This time you got ' + str(newpoints) + ' points; for ' + str(newtotal) + ' total.', 1, (50, 50, 50))
        self.screen.blit(txtFPS, (100,300))
        self.screenFlip()

    def drawDot(self, obj, color=(0,0,0)):
        cc = [[]]*3
        for i in range(len(color)):
            cc[i] = int(color[i]*255)
        pygame.draw.circle(self.screen, (0,0,0), (obj[0]*self.deg2pix+self.center[0], obj[1]*self.deg2pix+self.center[1]), 12)
        pygame.draw.circle(self.screen, cc, (obj[0]*self.deg2pix+self.center[0], obj[1]*self.deg2pix+self.center[1]), 10)

    def drawTarget(self, obj, color=(0,0,0)):
        cc = [[]]*3
        for i in range(len(color)):
            cc[i] = int(color[i]*255)
        pygame.draw.circle(self.screen, (0,0,0), (obj[0]*self.deg2pix+self.center[0], obj[1]*self.deg2pix+self.center[1]), 15)
        pygame.draw.circle(self.screen, (0,0,0), (obj[0]*self.deg2pix+self.center[0], obj[1]*self.deg2pix+self.center[1]), 12)
        pygame.draw.circle(self.screen, cc, (obj[0]*self.deg2pix+self.center[0], obj[1]*self.deg2pix+self.center[1]), 10)

    def drawState(self, states, Ts, centervoid=0):
        self.screenFill()
        if centervoid > 0:
            pygame.draw.circle(self.screen, (0,0,0), (self.center[0], self.center[1]), centervoid*self.deg2pix, 2)
        # self.writeInstructions(nt,no)
        self.drawFixation((0,0))
        states.reverse()
        for S in states:
            if S[4] not in Ts:
                col = S[5]
                self.drawDot(S, colorsys.hls_to_rgb(col, 0.5, 0.5))
        for S in states:
            if S[4] in Ts:
                col = S[5]
                self.drawTarget(S, colorsys.hls_to_rgb(col, 0.5, 0.5))
        self.screenFlip()

    def drawFixation(self, fixation, color=(0,0,0)):
        pygame.draw.circle(self.screen, color, (fixation[0]*self.deg2pix+self.center[0], fixation[1]*self.deg2pix+self.center[1]), 5)

    def screenFill(self, color=(255,255,255)):
        self.screen.fill(color)

    def screenFlip(self):
        pygame.display.flip()

    def close(self):
        pygame.display.quit()
