jQuery(document).ready(function() {
	var lastTF = '';
	for (i = 0; i < result.length; i++) {
		if (lastTF != result[i].TF) {
			lastTF = result[i].TF;
			var title = "<tr><td class='table-title' colspan='4' onclick='showConsensus(\""+result[i].Organism+"\", \""+result[i].TF+"\")'>"+result[i].Organism+"("+result[i].TF+")</td></tr>";
			jQuery("#result").append(title);
		}
		var text = "<tr class='"+(i%2==0?"even":"odd")+"'>"+
			"<td class='clicker', onclick='scrollTo(\""+result[i].Sequence+result[i].Position+result[i].Direction+"\");'>"+result[i].Sequence+"</td>"+
			"<td>"+result[i].Value+"</td>"+
			"<td>"+result[i].Position+"</td>"+
			"<td>"+result[i].Direction+"</td>"+                                                    
			   "</tr>";
		jQuery("#result").append(text);
	}
	
	for (j = 0; j < result.length; j++) {
		jQuery("#"+result[j].Sequence+result[j].Position+result[j].Direction).cluetip({splitTitle: '|', showTitle: false, cluetipClass: 'rounded', dropShadow: false, arrows: true});
	}

	jQuery("#dialog").dialog({autoOpen: false, resizable: false, width: 380, position: top, modal: true});
});

function scrollTo(id) {
	jQuery('html, body').animate({
		scrollTop: jQuery("#"+id+"").offset().top
	}, 2000);
}

function setModel() {
	var checkModels = jQuery("input:checked");
	var setModels = new Object();
	for(i=0; i < checkModels.length; i++) {
		if (setModels[checkModels[i].name] == null) {
			setModels[checkModels[i].name] = new Array();
		} 
		setModels[checkModels[i].name].push(checkModels[i].id);
	}

	var orgNames = jQuery("label[name=model]");
	for (i = 0; i < orgNames.length; i++) {
		var selected = setModels[orgNames[i].id];
		jQuery("#"+orgNames[i].id+"").html((selected != null ? '[' + setModels[orgNames[i].id].join() + ']' : ''));
	}
}

function showConsensus(org, tf) {
	jQuery("#dialog").dialog('option', 'title', org + '(' + tf+  ')'); 
	jQuery("#dialog").dialog("open");
	return false;
}
