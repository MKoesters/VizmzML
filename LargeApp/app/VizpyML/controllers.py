from flask import Blueprint, request, render_template, \
                  flash, g, session, redirect, url_for

Vizpymzml = Blueprint("Vizpymzml", __name__, url_prefix='/VizML')

@Vizpymzml.route('/')
def index():
	return render_template('VizpyML/index.html')

