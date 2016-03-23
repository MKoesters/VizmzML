from flask import Blueprint

Vizpymzml = Blueprint("Vizpymzml", __name__)

@main.route('/')
def index():
	return render_templat('')

