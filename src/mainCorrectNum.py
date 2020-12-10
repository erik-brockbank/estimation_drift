import pygame
import math, sys, os, shutil
from experimenter import *
from numpy import random as rand
import random
import time

white = (255,255,255)
black = (0,0,0)

diameter = 20
samples = 100
likelihood = .02
#calibrations = 3

def instruct(drawer, text):
    # display instructions onscreen, then wait for space or return to continue.
    drawer.screenFill((white))
    drawer.writeInstructions(text)
    drawer.screenFlip()
    while True:
        event = pygame.event.wait()
        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_ESCAPE:
                break
            if event.key == pygame.K_SPACE or event.key == pygame.K_RETURN:
                drawer.screen.fill(white)
                pygame.display.flip()
                return
    pygame.quit()
    sys.exit(0)

def geometricSample():
    ndots = []
    while len(ndots) < samples:
        ndot = rand.geometric(likelihood,1)[0]
        if ndot > 2:
            ndots.append(ndot)
    return ndots

def buildConfiguration(n):
    configuration = []
    for i in range(1,n):
        x = random.randint(int(.05*dispSize[0]), dispSize[0] - int(.05*dispSize[0]) - diameter)
        y = random.randint(int(.05*dispSize[1]), dispSize[1] - int(.05*dispSize[1]) - diameter)
        rect = pygame.Rect(x,y,diameter,diameter)
        if rect.collidelist(configuration) == -1:
            configuration.append(rect)
    return configuration

def buildConfigurations():
    configurations = []
    ndots = geometricSample()
    for n in ndots:
        #for calibration in range(0,calibrations):
        configuration = buildConfiguration(n)
        configurations.append(configuration)
    return configurations

def drawDots(drawer, configuration):
    # makes dots
    drawer.screenFill((white))
    drawer.screenFlip()
    pygame.time.wait(250)
    for rect in configuration:
        # pygame.draw.circle(drawer.screen, (255,0,0), rect.center, diameter/2, 0)
        pygame.draw.circle(drawer.screen, (255,0,0), rect.center, int(diameter/2), 0) #erikb
    drawer.screenFlip()
    pygame.time.wait(250)
    drawer.screenFill((white))
    drawer.screenFlip()

def runExperiment(f):
    feedback = False
    parts = 5
    a_choices = [.75,1.25]
    for part in range(0,parts):
        tot_points = 0
        configurations = buildConfigurations()
        if feedback == False:
            a = 0
            points = 0
            for configuration in configurations:
                drawDots(drawer,configuration)
                start_time = time.time()
                displayText ('Guess how many dots were on the screen', 200, 50, black)
                answer = inputValue()
                resp_time = time.time() - start_time
                f.write('%d,%d,%d,%r,%f,%d,%.2f\n' % (part+1,len(configuration),answer,feedback,a,points,resp_time))
        if feedback == True:
            a = random.choice(a_choices)
            a_choices.remove(a)
            for configuration in configurations:
                drawDots(drawer,configuration)
                start_time = time.time()
                start_time = time.time()
                displayText ('Guess how many dots were on the screen', 200, 50, black)
                answer = inputValue()
                resp_time = time.time() - start_time
                if len(configuration) <= 30:
                    points = max(0, round(10 - 10*(abs(math.log10(len(configuration)) - math.log10(answer)))))
                if len(configuration) > 30:
                    points = (round(10**(a*(math.log10(answer) - math.log10(30)) + math.log10(30))))/10
                print(points) #eb
                tot_points += points
                instruct(drawer, ['You scored %d points for a total of %d!' % (points, tot_points),
                                  'press ENTER or SPACE to continue'])
                f.write('%d,%d,%d,%r,%f,%d,%.2f\n' % (part+1,len(configuration),answer,feedback,a,points,resp_time))
        if part+1 == parts and feedback == False:
            instruct(drawer, ['You are done with the experiment!',
                              'press ENTER or SPACE to exit.'])
        if (part+1 != parts) and feedback == True:
            instruct(drawer, ['You are done with %d/%d of the experiment!' % (part+1,parts),
                              'In this section you have earned: %d points.' % (tot_points),
                              'You can now briefly stretch and then',
                              'press ENTER or SPACE to continue'])
        if (part+1 != parts) and feedback == False:
            instruct(drawer, ['You are done with %d/%d of the experiment!' % (part+1,parts),
                              'You can now briefly stretch and then',
                              'press ENTER or SPACE to continue'])
        feedback = not feedback

def checkInput(instr, nextchar):
    if len(instr) != 5:
        return True
    else:
        return False
    return True

def inputValue():
    # Set up the string
    instr = ''
    font = pygame.font.Font(None, 55)

    # Define the width of the text box:
    width = font.size('9999999')[0]

    # Draw a white background box for the rectangle
    bkrect = pygame.Rect(screen.get_rect().width/2 - width/5, screen.get_rect().height/2, width, font.get_linesize())
    screen.fill(white, bkrect)

    # Draw the text onto the rectangle
    txt = font.render('#' + instr, True, black)
    screen.blit(txt, bkrect)

    # Get alphanumeric key presses
    while True:
        for event in pygame.event.get():
            if event.type == pygame.KEYDOWN:
                nextchar = pygame.key.name(event.key)
                if event.key == event.key >= pygame.K_0 and event.key <= pygame.K_9:
                    if checkInput(instr, nextchar):
                        instr += nextchar
                elif event.key == pygame.K_BACKSPACE or event.key == pygame.K_DELETE:
                    if len(instr) > 0:
                        instr = instr[:-1]
                elif event.key == pygame.K_RETURN or event.key == pygame.K_SPACE:
                    if instr == '':
                        pass
                    else:
                        val = int(instr)
                        return val
                elif event.key == pygame.K_ESCAPE:
                    pygame.display.quit()
                    sys.exit(0)
            if event.type == pygame.QUIT:
                pygame.display.quit()
                sys.exit(0)

        # Redraw the text
        bkrect = pygame.Rect(screen.get_rect().width/2 - width/5, screen.get_rect().height/2, width, font.get_linesize())
        screen.fill(white, bkrect)

        # Draw the text onto the rectangle
        txt = font.render('# ' + instr, True, black)
        screen.blit(txt, bkrect)
        pygame.display.flip()


def displayText(text, top, size, color):
    drawer.screenFill((white))
    drawer.screenFlip()
    font = pygame.font.Font(None, size)
    # rectangular bitmap for text
    textRect = font.render(text, True, color)
    # where the text rectangle will be flipped
    textpos = textRect.get_rect(centerx = screen.get_rect().centerx)
    textpos.top = top
    # Block image transfer (combines 2 bitmaps : screen and textRect)
    screen.blit(textRect, textpos)
    pygame.display.flip()

# Initialize the game engine
pygame.init()
# Get subject ID and set up a file
# subject = raw_input('Enter subject ID here: ')
subject = input('Enter subject ID here: ') #eb
print('Subject ID is:', subject) # eb
# Set up output
f = open(subject + '.csv', 'w')
fname = subject + '.csv'
f.write('part,num_dots,answer,feedback,a,points,time\n')
# Set up output directory
cwd = os.getcwd()
dst = cwd + '/OUTPUT/'
# Initialize drawer
drawer = Experimenter()
# Set up the screen
screen=pygame.display.set_mode((0,0), pygame.FULLSCREEN)

dispSize = [int(drawer.screen.get_width()*.95), int(drawer.screen.get_height()*.95)]
# start timing the experiment
stime = time.time()
instruct(drawer, ['Welcome to the Experiment!',
                  '',
                  'You are going to be flashed some number of dots on the screen,',
                  'and your objective will be to make a guess about how many you saw.',
                  'There are going to be several sections and on some of them you will be awarded points.',
                  'The point system is based on the accuracy of both guesses so try your best.',
                  '',
                  'press ENTER or SPACE to start.'])

runExperiment(f)
# how long does the experiment take
print("expt done; time = %.2f minutes" % (int(time.time() - stime)/60.0)) # eb

pygame.quit()
sys.exit(0)
# Movie the file into the Output folder when the experiment is over.
shutil.move(cwd + fname , dst)
